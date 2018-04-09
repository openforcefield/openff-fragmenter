import os
import unittest
from pkg_resources import resource_filename

try:
    import openeye.oechem
    has_openeye = True
except ImportError:
    has_openeye = False


class FileIOTestCase(unittest.TestCase):

    def setUp(self):
        self._empty_writes()
        try:
            os.makedirs(get_fn('writes'))
        except OSError:
            pass

    def tearDown(self):
        self._empty_writes()
        try:
            os.rmdir(get_fn('writes'))
        except OSError:
            pass

    def _empty_writes(self):
        """ Empty the "writes" directory """
        try:
            for f in os.listdir(get_fn('writes')):
                os.unlink(get_fn(f, written=True))
        except OSError:
            pass

    def get_writes_dir(self):
        write_dir = resource_filename('torsionfit', os.path.join('tests', 'reference', 'writes'))
        if not os.path.exists(write_dir):
            self.setUp()
        return write_dir


def get_fn(filename, written=False):
    """Get the full path to one of the reference files shipped for testing

        These files are in torsionfit/testing/reference

    :param
        name: str
            Name of file to load

    :returns
        fn : str
            full path to file
    """
    if written:
        fn = resource_filename('torsionfit', os.path.join('tests', 'reference', 'writes', filename))
    else:
        fn = resource_filename('torsionfit', os.path.join('tests', 'reference', filename))

    #if not os.path.exists(fn):
    #    raise ValueError('%s does not exist. If you just added it you will have to re install' % fn)

    return fn



