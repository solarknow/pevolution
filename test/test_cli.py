import unittest

from helpers.cli import parse_bool_arg


class TestCli(unittest.TestCase):
    def test_parse_bool_arg_truthy_inputs(self):
        for value in ["1", "true", "TRUE", " yes ", "Y"]:
            self.assertTrue(parse_bool_arg(value))

    def test_parse_bool_arg_falsey_inputs(self):
        for value in ["0", "false", "no", "n", "maybe", ""]:
            self.assertFalse(parse_bool_arg(value))


if __name__ == "__main__":
    unittest.main()
