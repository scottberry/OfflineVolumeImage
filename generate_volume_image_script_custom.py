#! /usr/bin/env python
import sys
import logging
import itertools

import os
from os.path import basename

import gc3libs
from gc3libs.cmdline import SessionBasedScript
from gc3libs import Application
from gc3libs.quantity import GB
from gc3libs.workflow import SequentialTaskCollection, ParallelTaskCollection

from tmclient import TmClient


#logger = logging.getLogger(__name__)
#logging.basicConfig(format='%(asctime)s %(filename)s %(funcName)s %(levelname)s: %(message)s',
#    level=logging.DEBUG, datefmt='%I:%M:%S',stream=sys.stdout)


if __name__ == '__main__':
    from generate_volume_image_script_custom import GenerateVolumeImageScriptCustom
    GenerateVolumeImageScriptCustom().run()


class GenerateVolumeImageScriptCustom(SessionBasedScript):
    '''
    Script to generate volume images for a whole experiment using offline
    beads images
    '''

    def __init__(self):
        super(GenerateVolumeImageScriptCustom, self).__init__(version='1.0')

    def setup_args(self):
        self.add_param('--host', type=str, help=('TissueMAPS host address'))
        self.add_param('--username', type=str, help=('TissueMAPS username'))
        self.add_param('--password', type=str, help=('TissueMAPS password'))
        self.add_param('--plate', type=str, help=('Plate name'))
        self.add_param('--experiment', type=str, nargs='+',
                       help=('TissueMAPS experiment name'))
        self.add_param('--input_path', type=str, nargs='+', help=('Input path'))
        self.add_param('--output_path', type=str, nargs='+', help=('Output path'))
        self.add_param('--fname_stem', type=str, nargs='+', help=('Filename stem'))

    def new_tasks(self, extra):
        apps = [GenerateVolumeImagesSequentialPlates(self.params)]
        return apps

class GenerateVolumeImagesSequentialPlates(ParallelTaskCollection):
    '''
    Parallel task collection to loop over list of plates
    '''
    def __init__(self, params):
        self.params= params
        task_list = []

        for experiment, input_path, output_path, fname_stem in itertools.izip(params.experiment, params.input_path, params.output_path, params.fname_stem):
            print(experiment, input_path, output_path, fname_stem)
            self.params.experiment = experiment
            self.params.input_path = input_path
            self.params.output_path = output_path
            self.params.fname_stem = fname_stem
            task_list.append(
                GenerateVolumeImagesParallel(self.params)
            )

        ParallelTaskCollection.__init__(self, task_list, output_dir='')


class GenerateVolumeImagesParallel(ParallelTaskCollection):
    '''
    Run n_wells * n_sites instances of GenerateVolumeImagesApp in parallel
    '''

    def __init__(self, params):
        self.params = params
        task_list = []

        tmaps_api = TmClient(
            host=self.params.host,
            port=80,
            experiment_name=self.params.experiment,
            username=self.params.username,
            password=self.params.password
        )

        # find the site dimensions
        sites = tmaps_api.get_sites()
        for site in list(sites[1]):
            if site['plate_name'] == self.params.plate:
                well_name = site['well_name']
                x = site['x']
                y = site['y']

                task_list.append(
                    GenerateVolumeImagesApp(
                        self.params.host, self.params.username,
                        self.params.password, self.params.experiment,
                        self.params.plate, well_name, x, y,
                        self.params.input_path,
                        self.params.output_path,
                        self.params.fname_stem
                    )
                )
        ParallelTaskCollection.__init__(self, task_list, output_dir='')


class GenerateVolumeImagesApp(Application):
    '''
    Generate volume images for a single site
    '''

    def __init__(self, host, username, password, experiment,
                 plate, well_name, x, y, input_path,
                 output_path, fname_stem):

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
                '--plate_name', plate,
                '--well_name', well_name,
                '-x', x,
                '-y', y,
                '--input_path', input_path,
                '--output_path', output_path,
                '--fname_stem', fname_stem],
            inputs=['generate_volume_image_with_offline_beads.py'],
            outputs=[],
            output_dir=output_dir,
            stdout='stdout.txt',
            stderr='stderr.txt',
            requested_memory=3 * GB
        )
