"""
platform: win
env: chimera_env
name: chimera_run_ign.py
"""
import os
from config import *


def get_scripts(d):
    """
    get the script file(python src file)
    :param d: directory of scripts
    :return:
    """
    return [os.path.join(d, f) for f in os.listdir(d) if f[-2:] == 'py']


def main(sc):
    print('\n')
    scripts = get_scripts(sc)
    count = 0
    for script in scripts:
        count += 1
        os.system("chimera --nogui --silent --script {}".format(script))
        print(count)


script_dir = "D:\\AlexYu\\datasets\\dataset\\pdbbind2016\\example\\scripts"
main(script_dir)
