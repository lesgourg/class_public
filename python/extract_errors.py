# From the dumped stdout and stderr of a nosetests test_class.py, extract all
# the failed steps.
# Usage: python extract_errors.py output
from __future__ import print_function
import sys
import os


def main(path):
    """
    Create a shorter file containing only the errors from nosetests

    """
    assert os.path.isfile(path) is True

    trimmed_path = path + '_errors'
    destination = open(trimmed_path, 'w')
    contains_error = False
    with open(path, 'r') as source:
        text = source.readlines()
        start = 0
        for index, line in enumerate(text):
            if line.find('------------------') != -1:
                if text[index+2].find('----------------') != -1:
                    stop = index-1
                    # Check that an error is contained
                    if stop > 0:
                        for i in range(start, stop+1):
                            if text[i].startswith('E'):
                                contains_error = True
                    if contains_error:
                        print('Found an error')
                        for i in range(start, stop+1):
                            print(text[i], end=' ')
                            destination.write(text[i])
                    start = index
                    contains_error = False
                elif text[index+2].find('=================') != -1:
                    break
                else:
                    pass
    destination.close()


if __name__ == "__main__":
    print(sys.argv)
    if len(sys.argv) != 2:
        print('Please specify the output file to analyse')
        exit()
    else:
        main(sys.argv[-1])
