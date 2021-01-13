#!/usr/bin/python

"""
Translate an ASCII file in Oslo format from processSQL.py
to an acceptlist, for all files containing _oslo in the specified directory
"""

if __name__ == '__main__':

    import os,sys
    import scanf
    import secateur

    if len(sys.argv) != 2:
        print 'Usage: ./convertSQLtoAccept.py directory'
        sys.exit(1)

    dir=sys.argv[-1]
    assert(os.path.isdir(dir)), '%s is not a directory!' % (dir)
    secateur.ConvertChicagoFilesToAcceptlists(files=os.path.abspath(dir))

    sys.exit(0)
