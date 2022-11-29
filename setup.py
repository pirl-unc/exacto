# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


from os.path import dirname, join
from setuptools import setup, find_packages
from glob import glob


DIR = (dirname(__file__) or '.')


if __name__ == '__main__':
      setup(
            name='exacto',
            version='0.1.4.1',
            description='EXtracting And Counting Transcripts in Oncology.',
            author='Jin Seok (Andy) Lee',
            author_email='ajslee@unc.edu',
            install_requires=[
                  'pandas',
                  'pysam',
                  'pyensembl',
                  'pyvcf'
            ],
            packages=find_packages(),
            entry_points={
                  'console_scripts': [
                        'exacto=exactolib.cli.exacto_main:run'
                  ]
            }
      )
