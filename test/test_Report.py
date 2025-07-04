import unittest


class TestReport(unittest.TestCase):
    def setUp(self):
        self.test_name = 'test'

    def test_something(self):
        self.assertEqual(True, False)  # add assertion here


if __name__ == '__main__':
    unittest.main()
