#! /usr/bin/env python

import os
from gc3libs.cmdline import SessionBasedScript
from gc3libs import Application
from gc3libs.quantity import GB
from gc3libs.workflow import ParallelTaskCollection

from tmclient import TmClient


if __name__ == '__main__':
    from generate_volume_image_script import GenerateVolumeImageScript
    GenerateVolumeImageScript().run()


class GenerateVolumeImageScript(SessionBasedScript):
    '''
    Script to generate volume images for a whole experiment using offline
    beads images
    '''

    def __init__(self):
        super(GenerateVolumeImageScript, self).__init__(version='1.0')

    def setup_args(self):
        self.add_param(
            '-H', '--host', default='app.tissuemaps.org',
            help='name of TissueMAPS server host'
        )
        self.add_param(
            '-P', '--port', type=int, default=80,
            help='number of the port to which the server listens (default: 80)'
        )
        self.add_param(
            '-u', '--user', dest='username', required=True,
            help='name of TissueMAPS user'
        )
        self.add_param(
            '--password', required=True,
            help='password of TissueMAPS user'
        )
        self.add_param(
            '-e', '--experiment', required=True,
            help='experiment name'
        )
        self.add_param(
            '-p', '--plate_name', type=str, default='plate01',
            help='plate name'
        )
        self.add_param(
            '--threshold', type=int, default=12,
            help='threshold for bead detection'
        )
        self.add_param(
            '--mean_size', type=int, default=5,
            help='mean size of bead'
        )
        self.add_param(
            '--min_size', type=int, default=7,
            help='min size of bead'
        )
        self.add_param(
            '--filter_type', type=str, default='log_2d',
            help='min size of bead'
        )
        self.add_param(
            '--minimum_bead_intensity', type=int, default=135,
            help='minimum_bead_intensity'
        )
        self.add_param(
            '--z_step', type=float, default=0.4,
            help='z-step in microns'
        )
        self.add_param(
            '--pixel_size', type=float, default=0.1625,
            help='pixel size in microns'
        )
        self.add_param(
            '--alpha', type=int, default=250,
            help='alpha_shape smoothing parameter'
        )
        self.add_param(
            '--smooth', type=float, default=1.2,
            help='surface volume image smoothing parameter'
        )
        self.add_param(
            '--input_path', type=str, required=True,
            help='path to directory containing bead images'
        )
        self.add_param(
            '--output_path', type=str, required=True,
            help='path to write volume image'
        )
        self.add_param(
            '--channel_string', type=str,
            default='C03',
            help=('channel number for beads (e.g. C03)')
        )

    def new_tasks(self, extra):
        apps = [GenerateVolumeImagesParallel(self.params)]
        return apps


class GenerateVolumeImagesParallel(ParallelTaskCollection):
    '''
    Run n_wells * n_sites instances of GenerateVolumeImagesApp in parallel
    '''

    def __init__(self, params):
        self.params = params
        task_list = []

        print("GenerateVolumeImagesParallel: contacting tissumaps server")
        tmaps_api = TmClient(
            host=self.params.host,
            port=80,
            experiment_name=self.params.experiment,
            username=self.params.username,
            password=self.params.password
        )

        # find the site dimensions
        sites = tmaps_api.get_sites()
        print("found %d sites on tissuemaps",len(sites))
        for site in sites:
            if site['plate_name'] == self.params.plate_name:
                well_name = site['well_name']
                x = site['x']
                y = site['y']

                task_list.append(
                    GenerateVolumeImagesApp(
                        self.params.host,
                        self.params.username,
                        self.params.password,
                        self.params.experiment,
                        self.params.plate_name,
                        well_name,
                        x,
                        y,
                        self.params.threshold,
                        self.params.mean_size,
                        self.params.min_size,
                        self.params.filter_type,
                        self.params.minimum_bead_intensity,
                        self.params.z_step,
                        self.params.pixel_size,
                        self.params.alpha,
                        self.params.smooth,
                        self.params.input_path,
                        self.params.output_path,
                        self.params.channel_string
                    )
                )
        ParallelTaskCollection.__init__(self, task_list, output_dir='')


class GenerateVolumeImagesApp(Application):
    '''
    Generate volume images for a single site
    '''

    def __init__(self, host, username, password, experiment,
                 plate_name, well_name, x, y,
                 threshold, mean_size, min_size, filter_type,
                 minimum_bead_intensity,
                 z_step, pixel_size,
                 alpha, smooth,
                 input_path, output_path, channel_string):

        out = 'vol_img_{w}_x{x:02d}_y{y:02d}'.format(x=x,y=y,w=well_name)
        output_dir = os.path.join(experiment, out)
        Application.__init__(
            self,
            arguments=[
                './generate_volume_image_with_offline_beads.py',
                '-H', host,
                '-u', username,
                '--password', password,
                '-e', experiment,
                '--plate_name', plate_name,
                '--well_name', well_name,
                '-x', x,
                '-y', y,
                '--threshold', threshold,
                '--mean_size', mean_size,
                '--min_size', min_size,
                '--filter_type', filter_type,
                '--minimum_bead_intensity', minimum_bead_intensity,
                '--z_step', z_step,
                '--pixel_size', pixel_size,
                '--alpha', alpha,
                '--smooth', smooth,
                '--input_path', input_path,
                '--output_path', output_path,
                '--channel_string', channel_string],
            inputs=['generate_volume_image_with_offline_beads.py'],
            outputs=[],
            output_dir=output_dir,
            stdout='stdout.txt',
            stderr='stderr.txt',
            requested_memory=3 * GB
        )
