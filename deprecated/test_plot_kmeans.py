import unittest
import sys
import random
from plot_kmeans import import_file

class TestPlotKMeans(unittest.TestCase):

    def test_import_file(self):
        f = open('test_import.fa', 'a')
        for i in range(1, 100):
            f.write('>test{}'.format(counter))
            f.write('AABBCCDD')
        f.close()
        imported_proteins = import_file(f)
        r = imported_proteins[0][1]
        os.remove()
        s = 'AABBCCDD'
        self.assertEqual(r, s)

            
if __name__ == '__main__':
    unittest.main()