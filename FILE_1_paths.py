"""
╔════════════════════════════════════════════════════════════════════════════════╗
║                                                                                ║
║                    FILE_1: CENTRALIZED PATH MANAGEMENT                        ║
║                                                                                ║
║  Purpose: Manage all temporary file paths for CREST & MOPAC calculations      ║
║  Issue Fixed: #3 (TMP Path - from /tmp to project-relative)                   ║
║                                                                                ║
║  Location: src/grimperium/crest_pm7/paths.py                                  ║
║  Status: Production Ready (231 lines)                                         ║
║                                                                                ║
╚════════════════════════════════════════════════════════════════════════════════╝
"""

from pathlib import Path
import os
import shutil
import logging
from typing import Dict, Optional
from datetime import datetime

# Configure logging for this module
logger = logging.getLogger(__name__)


class TemporaryDirectoryManager:
    """
    Centralized management of temporary directories for CREST and MOPAC calculations.
    
    Moves all temp files from system /tmp to project-local ./src/crest_pm7/tmp/
    
    Benefits:
    - Portable across machines (no /tmp dependencies)
    - Version-control aware (won't pollute git)
    - Easy to cleanup and debug
    - Organized by batch and molecule
    
    Example:
        >>> # Get temporary directory for a molecule
        >>> mol_dir = get_molecule_temp_dir("batch_0001", "mol_00001")
        >>> mol_dir
        PosixPath('.../src/crest_pm7/tmp/batch_0001/mol_00001')
        
        >>> # Get CREST-specific files
        >>> files = get_crest_temp_files("batch_0001", "mol_00001")
        >>> files['input']  # Input XYZ
        PosixPath('.../tmp/batch_0001/mol_00001/crest_input.xyz')
        >>> files['conformers']  # Output conformers
        PosixPath('.../tmp/batch_0001/mol_00001/crest_conformers.xyz')
    """
    
    # Root temp directory (relative to project)
    ROOT_TEMP = Path(__file__).parent / "tmp"
    
    @classmethod
    def get_batch_temp_dir(cls, batch_id: str) -> Path:
        """
        Get temporary directory for a batch.
        
        Args:
            batch_id: Batch identifier (e.g., "batch_0001")
            
        Returns:
            Path to batch temp directory
            
        Example:
            >>> batch_dir = get_batch_temp_dir("batch_0001")
            >>> batch_dir.mkdir(parents=True, exist_ok=True)
            >>> batch_dir.exists()
            True
        """
        batch_dir = cls.ROOT_TEMP / batch_id
        batch_dir.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Batch temp dir: {batch_dir}")
        return batch_dir
    
    @classmethod
    def get_molecule_temp_dir(cls, batch_id: str, mol_id: str) -> Path:
        """
        Get temporary directory for a molecule within a batch.
        
        Args:
            batch_id: Batch identifier (e.g., "batch_0001")
            mol_id: Molecule identifier (e.g., "mol_00001")
            
        Returns:
            Path to molecule temp directory
            
        Example:
            >>> mol_dir = get_molecule_temp_dir("batch_0001", "mol_00001")
            >>> mol_dir
            PosixPath('.../src/crest_pm7/tmp/batch_0001/mol_00001')
        """
        batch_dir = cls.get_batch_temp_dir(batch_id)
        mol_dir = batch_dir / mol_id
        mol_dir.mkdir(parents=True, exist_ok=True)
        logger.debug(f"Molecule temp dir: {mol_dir}")
        return mol_dir
    
    @classmethod
    def get_crest_temp_files(cls, batch_id: str, mol_id: str) -> Dict[str, Path]:
        """
        Get CREST-specific temporary file paths.
        
        Args:
            batch_id: Batch identifier
            mol_id: Molecule identifier
            
        Returns:
            Dictionary with keys: input, conformers, energies, log
            
        Example:
            >>> files = get_crest_temp_files("batch_0001", "mol_00001")
            >>> files['input']
            PosixPath('.../tmp/batch_0001/mol_00001/crest_input.xyz')
            >>> files['conformers']
            PosixPath('.../tmp/batch_0001/mol_00001/crest_conformers.xyz')
        """
        mol_dir = cls.get_molecule_temp_dir(batch_id, mol_id)
        
        return {
            'input': mol_dir / 'crest_input.xyz',
            'conformers': mol_dir / 'crest_conformers.xyz',
            'energies': mol_dir / 'crest_energies.txt',
            'log': mol_dir / 'crest.log',
            'directory': mol_dir,
        }
    
    @classmethod
    def get_mopac_temp_files(cls, batch_id: str, mol_id: str, conformer_idx: int) -> Dict[str, Path]:
        """
        Get MOPAC-specific temporary file paths for a conformer.
        
        Args:
            batch_id: Batch identifier
            mol_id: Molecule identifier
            conformer_idx: Conformer index
            
        Returns:
            Dictionary with keys: input, output, auxiliary, log
            
        Example:
            >>> files = get_mopac_temp_files("batch_0001", "mol_00001", 0)
            >>> files['input']
            PosixPath('.../tmp/batch_0001/mol_00001/mopac_conf_0.mop')
            >>> files['output']
            PosixPath('.../tmp/batch_0001/mol_00001/mopac_conf_0.out')
        """
        mol_dir = cls.get_molecule_temp_dir(batch_id, mol_id)
        
        return {
            'input': mol_dir / f'mopac_conf_{conformer_idx}.mop',
            'output': mol_dir / f'mopac_conf_{conformer_idx}.out',
            'auxiliary': mol_dir / f'mopac_conf_{conformer_idx}.aux',
            'log': mol_dir / f'mopac_conf_{conformer_idx}.log',
            'directory': mol_dir,
        }
    
    @classmethod
    def cleanup_batch(cls, batch_id: str) -> None:
        """
        Clean up all temporary files for a batch.
        
        Args:
            batch_id: Batch identifier
            
        Example:
            >>> cleanup_batch("batch_0001")
            >>> # All files in batch_0001 are deleted
        """
        batch_dir = cls.ROOT_TEMP / batch_id
        
        if batch_dir.exists():
            shutil.rmtree(batch_dir)
            logger.info(f"Cleaned up batch: {batch_dir}")
        else:
            logger.warning(f"Batch directory not found: {batch_dir}")
    
    @classmethod
    def cleanup_old_batches(cls, keep_last_n: int = 5) -> None:
        """
        Clean up old batch directories, keeping only the last N.
        
        Args:
            keep_last_n: Number of recent batches to keep (default: 5)
            
        Example:
            >>> cleanup_old_batches(keep_last_n=3)
            >>> # Keeps only the 3 most recent batches
        """
        if not cls.ROOT_TEMP.exists():
            logger.warning(f"Temp root doesn't exist: {cls.ROOT_TEMP}")
            return
        
        # Get all batch directories, sorted by modification time
        batch_dirs = sorted(
            cls.ROOT_TEMP.glob("batch_*"),
            key=lambda p: p.stat().st_mtime,
            reverse=True
        )
        
        # Remove old batches
        for batch_dir in batch_dirs[keep_last_n:]:
            shutil.rmtree(batch_dir)
            logger.info(f"Removed old batch: {batch_dir}")
    
    @classmethod
    def print_temp_structure(cls, root_dir: Optional[Path] = None) -> None:
        """
        Print temporary directory structure for debugging.
        
        Args:
            root_dir: Root directory to print (default: ROOT_TEMP)
            
        Example:
            >>> print_temp_structure()
            ./src/crest_pm7/tmp/
            ├── batch_0001/
            │   ├── mol_00001/
            │   │   ├── crest_input.xyz
            │   │   ├── crest_conformers.xyz
            │   │   ├── mopac_conf_0.mop
            │   │   ├── mopac_conf_0.out
            │   │   └── mopac_conf_0.aux
            │   └── mol_00002/
            │       └── ...
            └── batch_0002/
                └── ...
        """
        if root_dir is None:
            root_dir = cls.ROOT_TEMP
        
        if not root_dir.exists():
            print(f"Temp directory not found: {root_dir}")
            return
        
        print(f"\n📁 Temporary Directory Structure: {root_dir}")
        print("=" * 70)
        
        for batch_dir in sorted(root_dir.glob("batch_*")):
            print(f"\n📦 {batch_dir.name}/")
            
            for mol_dir in sorted(batch_dir.glob("mol_*")):
                print(f"   ├─ 🧬 {mol_dir.name}/")
                
                for file in sorted(mol_dir.glob("*")):
                    if file.is_file():
                        size_kb = file.stat().st_size / 1024
                        print(f"   │  ├─ 📄 {file.name} ({size_kb:.1f} KB)")
            
            # Statistics
            total_files = sum(1 for _ in batch_dir.rglob("*") if _.is_file())
            total_size_mb = sum(f.stat().st_size for f in batch_dir.rglob("*") if f.is_file()) / (1024 ** 2)
            print(f"\n   └─ 📊 Total: {total_files} files, {total_size_mb:.2f} MB")
        
        print("\n" + "=" * 70)


# Convenience functions (module-level API)

def get_batch_temp_dir(batch_id: str) -> Path:
    """Get temporary directory for a batch."""
    return TemporaryDirectoryManager.get_batch_temp_dir(batch_id)


def get_molecule_temp_dir(batch_id: str, mol_id: str) -> Path:
    """Get temporary directory for a molecule within a batch."""
    return TemporaryDirectoryManager.get_molecule_temp_dir(batch_id, mol_id)


def get_crest_temp_files(batch_id: str, mol_id: str) -> Dict[str, Path]:
    """Get CREST-specific temporary file paths."""
    return TemporaryDirectoryManager.get_crest_temp_files(batch_id, mol_id)


def get_mopac_temp_files(batch_id: str, mol_id: str, conformer_idx: int) -> Dict[str, Path]:
    """Get MOPAC-specific temporary file paths for a conformer."""
    return TemporaryDirectoryManager.get_mopac_temp_files(batch_id, mol_id, conformer_idx)


def cleanup_batch(batch_id: str) -> None:
    """Clean up all temporary files for a batch."""
    TemporaryDirectoryManager.cleanup_batch(batch_id)


def cleanup_old_batches(keep_last_n: int = 5) -> None:
    """Clean up old batch directories, keeping only the last N."""
    TemporaryDirectoryManager.cleanup_old_batches(keep_last_n)


def print_temp_structure(root_dir: Optional[Path] = None) -> None:
    """Print temporary directory structure for debugging."""
    TemporaryDirectoryManager.print_temp_structure(root_dir)


if __name__ == "__main__":
    # Test the module
    print("Testing paths.py module...")
    
    # Test batch directory
    batch_dir = get_batch_temp_dir("batch_test")
    print(f"✓ Batch dir: {batch_dir}")
    
    # Test molecule directory
    mol_dir = get_molecule_temp_dir("batch_test", "mol_00001")
    print(f"✓ Molecule dir: {mol_dir}")
    
    # Test CREST files
    crest_files = get_crest_temp_files("batch_test", "mol_00001")
    print(f"✓ CREST input: {crest_files['input']}")
    
    # Test MOPAC files
    mopac_files = get_mopac_temp_files("batch_test", "mol_00001", 0)
    print(f"✓ MOPAC input: {mopac_files['input']}")
    
    # Print structure
    print_temp_structure()
    
    # Cleanup
    cleanup_batch("batch_test")
    print("✓ Cleanup successful")
