import os
import sys

_PATH = os.path.expanduser('~') + '/fempkg/results/'


def make_dir(_name):
  dir = _name + '-0'
  i = 1
  while os.path.isdir(_PATH + dir):
    dir = dir.split('-')[0] + '-'+str(i)
    i+=1

  results_path = _PATH + dir
  os.mkdir(results_path)
  os.mkdir(results_path + '/vtk')
  return results_path

def sim_info_files(_path, _info):

  from datetime import datetime
  date = datetime.now().strftime("%H:%M:%S %d/%m/%Y")
  _file1 = open('/'.join((_path,'info.txt')), 'w+')
  _file1.write(' Simulation Parameters Info:\n\n')
  _file1.write(' Time started: {}\n\n'.format(date))
  for data in _info.items():
    _file1.write('{}: {}\n'.format(data[0], data[1]))
  _file1.close()

  from shutil import copy
  copy(sys.argv[0], _path)
  copy(_info['Mesh File'], _path)



if __name__ == '__main__':
  name = 'test'
  make_dir(name)
