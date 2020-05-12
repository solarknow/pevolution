import multiprocessing
import subprocess

jobs = open('jobs.txt').readlines()


def launch(query, label, domain):
    """
    Wrapper method to Main.py execution.
    """
    subprocess.call(['python', 'Main.py', query, label, domain])


if __name__ == '__main__':
    for i in jobs:
        spl = i.split()
        (name, quer) = (spl[0], spl[1])
        p = multiprocessing.Process(name=name, target=launch,
                                    args=(quer, name, 'all'))
        p.start()
