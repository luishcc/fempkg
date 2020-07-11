import os
import sys

_PATH = os.path.expanduser('~') + '/fempkg/results/'


def check_dir(_path):

    return os.path.isdir(_path)


def make_dir(_name):
    dir = _name + '-0'
    i = 1
    while check_dir(_PATH + dir):
        dir = dir.split('-')[0] + '-'+str(i)
        print(dir)
        i+=1
    os.mkdir(_PATH + dir)
    os.mkdir(_PATH + dir + '/vtk')


if __name__ == '__main__':
    name = 'test'
    make_dir(name)
