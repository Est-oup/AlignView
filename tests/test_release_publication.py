"""Release publication must not depend on TMPDIR sharing the output filesystem."""
import sys
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import build_release


class ReleasePublicationTests(unittest.TestCase):
    def test_replace_uses_a_stage_on_the_destination_filesystem(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);source=root/'scratch.whl';source.write_bytes(b'new artifact')
            dist=root/'dist';dist.mkdir();target=dist/'release.whl';target.write_bytes(b'old artifact')
            replace=build_release.os.replace
            def same_mount(staged,destination):
                self.assertTrue(staged.is_relative_to(destination.parent))
                return replace(staged,destination)
            with patch.object(build_release.os,'replace',side_effect=same_mount):
                build_release.publish_artifact(source,target)
            self.assertEqual(source.read_bytes(),b'new artifact')
            self.assertEqual(target.read_bytes(),b'new artifact')
            self.assertEqual(list(dist.iterdir()),[target])

    def test_failed_copy_leaves_previous_release_intact(self):
        with tempfile.TemporaryDirectory() as temporary:
            root=Path(temporary);source=root/'source';source.write_bytes(b'new')
            dist=root/'dist';dist.mkdir();target=dist/'release';target.write_bytes(b'old')
            def incomplete(source,staged):
                staged.write_bytes(b'incomplete')
                raise OSError('scratch storage failure')
            with patch.object(build_release.shutil,'copyfile',side_effect=incomplete):
                with self.assertRaises(OSError):build_release.publish_artifact(source,target)
            self.assertEqual(target.read_bytes(),b'old')
            self.assertEqual(list(dist.iterdir()),[target])


if __name__=='__main__':unittest.main()
