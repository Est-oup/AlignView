"""Preparation checks prevent stale or fictitious Bioconda source recipes."""
from __future__ import annotations
import hashlib
import io
import json
from pathlib import Path
import shutil
import sys
import tarfile
import tempfile
import unittest
from unittest.mock import patch

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
from prepare_bioconda import prepare,render_recipe,source_inventory,validate_public_url,verify_archive
from av_version import VERSION


class BiocondaPreparationTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder=Path(self.temp.name);self.source=self.folder/'source';self.source.mkdir()
        (self.source/'docs').mkdir()
        (self.source/'packaging/bioconda').mkdir(parents=True)
        (self.source/'alignview.py').write_text('def main(): return 0\n')
        (self.source/'av_version.py').write_text('VERSION = '+repr(VERSION)+'\n')
        (self.source/'pyproject.toml').write_text('[tool.setuptools]\npy-modules = ["alignview", "av_version"]\n')
        (self.source/'docs/README.md').write_text('Manual\n')
        (self.source/'LICENSE').write_text('MIT\n')
        (self.source/'README.md').write_text('AlignView\n')
        shutil.copy2(ROOT/'packaging/bioconda/meta.yaml.template',self.source/'packaging/bioconda/meta.yaml.template')
        self.archive=self.folder/'alignview.tar.gz'
        self.pack()
        self.sha=hashlib.sha256(self.archive.read_bytes()).hexdigest()

    def pack(self):
        with tarfile.open(self.archive,'w:gz') as archive:
            for path in source_inventory(self.source):
                archive.add(path,arcname='alignview-'+VERSION+'/'+path.relative_to(self.source).as_posix())

    def output(self,name='local'):
        return self.folder/name/'recipes/alignview'

    def test_local_candidate_has_real_checksum_entrypoint_and_installed_smoke_tests(self):
        output,evidence=prepare(archive=self.archive,local=True,maintainer='Est-oup',output=self.output(),root=self.source)
        source=(output/'meta.yaml').read_text()
        self.assertIn(self.sha,source)
        self.assertIn(self.archive.as_uri(),source)
        self.assertIn('alignview = alignview:main',source)
        self.assertIn('noarch: python',source)
        self.assertIn('pin_subpackage(name, max_pin="x")',source)
        self.assertIn('alignview --docs quickstart',source)
        self.assertIn('alignview report "bioconda smoke project"',source)
        self.assertNotIn('@@',source)
        data=json.loads(evidence.read_text())
        self.assertEqual(data['mode'],'local-candidate')
        self.assertFalse(data['conda_build_run'])
        self.assertIn('alignview.py',data['verified_files'])

    def test_stale_source_rejected_before_recipe_is_written(self):
        (self.source/'alignview.py').write_text('def main(): return 1\n')
        with self.assertRaisesRegex(ValueError,'differs from current sources'):
            prepare(archive=self.archive,local=True,maintainer='Est-oup',output=self.output(),root=self.source)
        self.assertFalse(self.output().exists())

    def test_public_recipe_download_must_match_independent_checksum_and_sources(self):
        def fetch(url,target):shutil.copyfile(self.archive,target)
        url='https://github.com/Est-oup/AlignView/releases/download/v'+VERSION+'/alignview-'+VERSION+'.tar.gz'
        with patch('prepare_bioconda.download',side_effect=fetch) as download:
            output,evidence=prepare(source_url=url,expected_sha256=self.sha,maintainer='Est-oup',output=self.output('public'),root=self.source)
        download.assert_called_once()
        self.assertIn(url,(output/'meta.yaml').read_text())
        self.assertEqual(json.loads(evidence.read_text())['mode'],'public-source-verified')
        with patch('prepare_bioconda.download',side_effect=fetch),self.assertRaisesRegex(ValueError,'SHA256'):
            prepare(source_url=url,expected_sha256='0'*64,maintainer='Est-oup',output=self.output('wrong'),root=self.source)
        self.assertFalse(self.output('wrong').exists())

    def test_unverifiable_submission_inputs_are_rejected(self):
        for url in ('http://example.org/alignview-'+VERSION+'.tar.gz',
                    'https://token@example.org/alignview-'+VERSION+'.tar.gz',
                    'https://example.org/main.tar.gz','https://example.org/refs/heads/'+VERSION+'.tar.gz',
                    'https://example.org/'+VERSION+'.tar.gz?token=abc'):
            with self.subTest(url=url),self.assertRaises(ValueError):validate_public_url(url)
        with self.assertRaisesRegex(ValueError,'requires --source-url'):
            prepare(archive=self.archive,maintainer='Est-oup',output=self.output(),root=self.source)
        with self.assertRaises(ValueError):render_recipe('url',self.sha,'not a handle',self.source)

    def test_traversal_links_and_duplicate_archive_members_are_rejected(self):
        for kind in ('traversal','link','duplicate'):
            with self.subTest(kind=kind):
                with tarfile.open(self.archive,'w:gz') as archive:
                    entry=tarfile.TarInfo('../outside' if kind=='traversal' else 'root/file')
                    if kind=='link':entry.type=tarfile.SYMTYPE;entry.linkname='/etc/passwd'
                    archive.addfile(entry,io.BytesIO(b''))
                    if kind=='duplicate':archive.addfile(entry,io.BytesIO(b''))
                with self.assertRaises(ValueError):verify_archive(self.archive,self.source)

    def test_existing_recipe_is_never_overwritten(self):
        output=self.output();output.mkdir(parents=True);(output/'keep').write_text('user file')
        with self.assertRaisesRegex(ValueError,'already exists'):
            prepare(archive=self.archive,local=True,maintainer='Est-oup',output=output,root=self.source)
        self.assertEqual((output/'keep').read_text(),'user file')


if __name__=='__main__':unittest.main()
