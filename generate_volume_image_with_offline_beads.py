#! /usr/bin/env python

import os
import re
import sys
import logging
import argparse
import numpy as np
import tifffile as tiff
import imageio
from tmclient import TmClient
from jtmodules import generate_volume_image as gvi
from jtmodules import project

logger = logging.getLogger(__name__)
logging.basicConfig(format='%(asctime)s %(filename)s %(funcName)s %(levelname)s: %(message)s',
    level=logging.DEBUG, datefmt='%I:%M:%S',stream=sys.stdout)


def parse_arguments():
    parser = argparse.ArgumentParser(
        prog='generate_volume_image_with_offline_beads',
        description=('Accesses DAPI and SE images from TissueMAPS instance'
                     ' to generate segmentation and generates a volume image'
                     ' from locally stored bead-stacks.')
    )
    parser.add_argument(
        '-v', '--verbosity', action='count', default=0,
        help='increase logging verbosity'
    )
    parser.add_argument(
        '-H', '--host', default='app.tissuemaps.org',
        help='name of TissueMAPS server host'
    )
    parser.add_argument(
        '-P', '--port', type=int, default=80,
        help='number of the port to which the server listens (default: 80)'
    )
    parser.add_argument(
        '-u', '--user', dest='username', required=True,
        help='name of TissueMAPS user'
    )
    parser.add_argument(
        '--password', required=True,
        help='password of TissueMAPS user'
    )
    parser.add_argument(
        '-e', '--experiment', required=True,
        help='experiment name'
    )
    parser.add_argument(
        '-p', '--plate_name', type=str, default='plate01',
        help='plate name'
    )
    parser.add_argument(
        '-w', '--well_name', type=str, required=True,
        help='well_name'
    )
    parser.add_argument(
        '-x', '--well_pos_x', type=int, required=True,
        help='well_pos_x'
    )
    parser.add_argument(
        '-y', '--well_pos_y', type=int, required=True,
        help='well_pos_y'
    )
    parser.add_argument(
        '--threshold', type=int, default=12,
        help='threshold for bead detection'
    )
    parser.add_argument(
        '--mean_size', type=int, default=5,
        help='mean size of bead'
    )
    parser.add_argument(
        '--min_size', type=int, default=7,
        help='min size of bead'
    )
    parser.add_argument(
        '--filter_type', type=str, default='log_2d',
        help='min size of bead'
    )
    parser.add_argument(
        '--minimum_bead_intensity', type=int, default=135,
        help='minimum_bead_intensity'
    )
    parser.add_argument(
        '--z_step', type=float, default=0.4,
        help='z-step in microns'
    )
    parser.add_argument(
        '--pixel_size', type=float, default=0.1625,
        help='pixel size in microns'
    )
    parser.add_argument(
        '--alpha', type=int, default=150,
        help='alpha_shape smoothing parameter'
    )
    parser.add_argument(
        '--smooth', type=int, default=1.2,
        help='surface volume image smoothing parameter'
    )
    parser.add_argument(
        '--input_path', type=str, required=True,
        help='path to directory containing bead images'
    )
    parser.add_argument(
        '--output_path', type=str, required=True,
        help='path to write volume image'
    )
    parser.add_argument(
        '--channel_string', type=str,
        default='C03',
        help=('channel number for beads (e.g. C03)')
    )

    return(parser.parse_args())


def main(args):

    tmaps_api = TmClient(
        host=args.host,
        port=args.port,
        experiment_name=args.experiment,
        username=args.username,
        password=args.password
    )

    dapi2D = tmaps_api.download_channel_image(
        channel_name='DAPI',
        plate_name=args.plate_name,
        well_name=args.well_name,
        well_pos_y=args.well_pos_y,
        well_pos_x=args.well_pos_x,
        cycle_index=0,
        tpoint=0,
        zplane=0,
        correct=True
    )

    se2D = tmaps_api.download_channel_image(
        channel_name='SE',
        plate_name=args.plate_name,
        well_name=args.well_name,
        well_pos_y=args.well_pos_y,
        well_pos_x=args.well_pos_x,
        cycle_index=0,
        tpoint=0,
        zplane=0,
        correct=True)

    # segment cells
    nuclei = _segmentPrimary(dapi2D)
    cells = _segmentSecondary(nuclei,se2D)

    wells = tmaps_api.get_wells()
    sites = tmaps_api.get_sites()

    # find the site dimensions
    for well in wells:
        if well['name'] == args.well_name:
            dimensions = well['dimensions']

    # convert the position into the Yokogawa field string
    field_str = str(args.well_pos_y * dimensions[1] + args.well_pos_x + 1).zfill(3)

    # find corresponding image files
    beads_files_all = [os.path.basename(full_path) for full_path in glob.glob(args.input_path + '*' + args.channel_string + '.png')]
    regex = (
        r'[^_]+_' + args.well_name +
        r'_T(?P<t>\d+)F' + field_str +
        r'L\d+A\d+Z(?P<z>\d+)' + args.channel_string +
        r'\.'
    )
    search_exp = re.compile(regex)
    beads_files_site = filter(search_exp.match, beads_files_all)

    for beads_file in beads_files_site:
        # allocate the data structure for beads
        beads3D = np.zeros(
            (sites[0]['height'], sites[0]['width'], len(beads_files_site)),
            dtype=np.uint16
        )
        logger.debug('load {0}'.format(beads_file))

        # get z-position and extension from filename
        matches = re.match(search_exp, beads_file)
        index = int(matches.group('z')) - 1
        file_type = os.path.splitext(beads_file)[1]

        # load image into array
        if file_type in set('tif','tiff','TIF','TIFF'):
            beads3D[:,:,index] = tiff.imread(
                os.path.join(args.input_path,beads_files_site)
            )
        elif file_type in set('.png','.PNG'):
            beads3D[:,:,index] = imageio.imread(
                os.path.join(args.input_path,beads_files_site)
            )

    beads3D = np.ascontiguousarray(beads3D)

    logger.debug('calculate volume image field {0}'.format(field_str))
    volume_image = gvi.main(
        image=beads3D,
        mask=cells.secondary_label_image,
        threshold=args.threshold,
        mean_size=args.mean_size,
        min_size=args.min_size,
        filter_type=args.filter_type,
        minimum_bead_intensity=args.minimum_bead_intensity,
        z_step=args.z_step,
        pixel_size=args.pixel_size,
        alpha=args.alpha,
        smooth=args.smooth,
        plot=False
    )

    vol_name = (args.experiment + '_' + args.well_name + '_T0001' +
        'F' + field_str + 'L01A02Z01C04.png')
    surface_name = (args.experiment + '_' + args.well_name + '_T0001' +
        'F' + field_str + 'L01A02Z01C05.png')
    max_proj_name = (args.experiment + '_' + args.well_name + '_T0001' +
        'F' + field_str + 'L01A02Z01C06.png')

    if not os.path.exists(args.output_path):
        os.makedirs(args.output_path)

    logger.debug('saving {0}'.format(vol_name))
    imageio.imsave(
        file=os.path.join(args.output_path,vol_name),
        data=volume_image.volume_image.astype(np.uint16)
    )

    logger.debug('saving {0}'.format(surface_name))
    imageio.imsave(
        file=os.path.join(args.output_path,surface_name),
        data=volume_image.smoothed_surface_image.astype(np.uint16)
    )

    logger.debug('making maximum projection field {0}'.format(field_str))
    max_proj = project.main(
        image=beads3D,
        method='max',
        plot=False
    )
    logger.debug('saving {0}'.format(max_proj_name))
    imageio.imsave(
        file=os.path.join(args.output_path,max_proj_name),
        data=max_proj.projected_image.astype(np.uint16)
    )

    return

def _segmentPrimary(dapi):
    from jtmodules import smooth, threshold_adaptive, fill, filter, label, register_objects, separate_clumps
    dapiSmooth = smooth.main(dapi, 'bilateral', 3)
    nucleiMask = threshold_adaptive.main(
        image=dapiSmooth.smoothed_image,
        method='niblack',
        kernel_size=211,
        constant=5,
        min_threshold=118,
        max_threshold=125
    )
    nucleiFilledMask = fill.main(nucleiMask.mask)
    nucleiFilledMaskSeparated = separate_clumps.main(
        mask=nucleiFilledMask.filled_mask,
        intensity_image=dapiSmooth.smoothed_image,
        min_area=5000,
        max_area=50000,
        min_cut_area=2000,
        max_circularity=0.7,
        max_convexity=0.92
    )
    nucleiFilledMaskFiltered = filter.main(
        mask=nucleiFilledMaskSeparated.separated_mask,
        feature='area',
        lower_threshold=1000,
        upper_threshold=None
    )
    nucleiLabel = label.main(mask=nucleiFilledMaskFiltered.filtered_mask)
    nuclei = register_objects.main(nucleiLabel.label_image)
    return nuclei


def _segmentSecondary(nuclei, celltrace):
    from jtmodules import segment_secondary, smooth
    celltraceSmooth = smooth.main(celltrace, 'bilateral', 5)
    cells = segment_secondary.main(nuclei.objects,
                                   celltraceSmooth.smoothed_image,
                                   contrast_threshold=5,
                                   min_threshold=118,
                                   max_threshold=125)
    return cells


if __name__ == '__main__':
    arguments = parse_arguments()
    main(arguments)
