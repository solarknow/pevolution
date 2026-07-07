import unittest
from unittest.mock import patch

from helpers import commands


class TestCommands(unittest.TestCase):
    @patch("helpers.commands.subprocess.run")
    def test_format_run_uses_check_true(self, mock_run):
        commands.format_run(["echo", "{value}"], value="ok")
        mock_run.assert_called_once_with(["echo", "ok"], check=True)

    @patch("helpers.commands.format_run")
    def test_run_phyml_uses_per_call_copy_without_mutating_global(self, mock_format_run):
        original_phyml = list(commands.PHYML)

        commands.run_phyml("infile", "JTT", "outfile", v="0.1")
        commands.run_phyml("infile", "JTT", "outfile", a="0.2")

        self.assertEqual(commands.PHYML, original_phyml)

        first_cmd = mock_format_run.call_args_list[0].args[0]
        second_cmd = mock_format_run.call_args_list[1].args[0]

        self.assertIn("-v", first_cmd)
        self.assertNotIn("-a", first_cmd)

        self.assertIn("-a", second_cmd)
        self.assertNotIn("-v", second_cmd)

    @patch("helpers.commands.format_run")
    def test_run_blast_builds_per_call_command_without_mutating_global(self, mock_format_run):
        original_blast = list(commands.BLAST)

        commands.run_blast("ABC123", 1e-5, "seed_org", "Homo sapiens", db="nr")
        commands.run_blast("ABC123", 1e-5, "seed_org", "Homo sapiens", db="local_db")

        self.assertEqual(commands.BLAST, original_blast)

        nr_cmd = mock_format_run.call_args_list[0].args[0]
        local_cmd = mock_format_run.call_args_list[1].args[0]

        self.assertIn("-remote", nr_cmd)
        self.assertIn("-entrez_query", nr_cmd)
        self.assertIn('"Homo sapiens[ORGN]"', nr_cmd)

        self.assertNotIn("-remote", local_cmd)
        self.assertNotIn("-entrez_query", local_cmd)


if __name__ == "__main__":
    unittest.main()
