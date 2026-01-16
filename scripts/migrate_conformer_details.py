"""Migration script to rename conformer_details JSON files from numeric to mol_XXXXX format."""

from pathlib import Path


def migrate_conformer_details_filenames():
    """Rename JSON files from numeric (2.json) to mol_XXXXX format (mol_00002.json)."""
    conformer_dir = Path("data/molecules_pm7/conformer_details")

    if not conformer_dir.exists():
        print(f"Directory not found: {conformer_dir}")
        return

    json_files = sorted(conformer_dir.glob("*.json"))
    migrated_count = 0

    for json_file in json_files:
        try:
            mol_id_num = int(json_file.stem)
        except ValueError:
            continue  # Already in new format or not numeric

        new_name = f"mol_{mol_id_num:05d}.json"
        new_path = conformer_dir / new_name

        json_file.rename(new_path)
        print(f"  Renamed: {json_file.name} -> {new_name}")
        migrated_count += 1

    print(f"\nMigrated {migrated_count} files")


if __name__ == "__main__":
    migrate_conformer_details_filenames()
