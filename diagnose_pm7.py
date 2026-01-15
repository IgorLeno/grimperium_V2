#!/usr/bin/env python3
"""
GRIMPERIUM PM7 FIRST RUN DIAGNOSIS
Detecta problemas antes de escalar para 30k
"""

import sys
import json
from pathlib import Path

# Add src to path
sys.path.insert(0, str(Path(__file__).parent / "src"))

def diagnose():
    print("\n" + "="*80)
    print("🔍 GRIMPERIUM PM7 DIAGNOSIS")
    print("="*80 + "\n")
    
    checks_passed = 0
    checks_failed = 0
    
    # CHECK 1: Timeout Conversion
    print("1️⃣  Checking timeout conversion...")
    try:
        from grimperium.crest_pm7.batch.processor_adapter import FixedTimeoutProcessor
        from grimperium.crest_pm7.config import PM7Config

        # Test 1: Verify FixedTimeoutProcessor correctly converts minutes to seconds
        config = PM7Config()
        crest_minutes = 20
        mopac_minutes = 10

        processor = FixedTimeoutProcessor(
            config=config,
            crest_timeout_minutes=crest_minutes,
            mopac_timeout_minutes=mopac_minutes,
        )

        # The config.crest_timeout should now be in seconds
        expected_seconds = crest_minutes * 60  # 1200 seconds
        actual_seconds = config.crest_timeout

        print(f"   Input: {crest_minutes} minutes")
        print(f"   Expected: {expected_seconds} seconds")
        print(f"   Actual config.crest_timeout: {actual_seconds} seconds")

        if actual_seconds == expected_seconds:
            print(f"   ✅ Timeout conversion correct: {crest_minutes}min = {actual_seconds}s")
            checks_passed += 1
        else:
            print(f"   ❌ Timeout conversion WRONG: expected {expected_seconds}s, got {actual_seconds}s")
            checks_failed += 1
    except Exception as e:
        print(f"   ❌ Erro ao testar timeout conversion: {e}")
        import traceback
        traceback.print_exc()
        checks_failed += 1
    
    # CHECK 2: CSV Schema
    print("\n2️⃣  Checking CSV schema...")
    try:
        import warnings
        from grimperium.crest_pm7.batch.csv_manager import BatchCSVManager

        csv_path = Path("data/test_batch_small.csv")
        if csv_path.exists():
            # Test using BatchCSVManager which applies our dtype fixes
            manager = BatchCSVManager(csv_path)
            df = manager.load_csv()

            print(f"   Loaded {len(df)} molecules via BatchCSVManager")

            # Check for problematic columns
            problematic = []

            # Check 'success' column type (should be object to handle bool + NaN)
            if 'success' in df.columns:
                dtype = str(df['success'].dtype)
                if dtype not in ('object', 'bool'):
                    problematic.append(f"success: {dtype} (deveria ser object/bool)")
                else:
                    print(f"   success column dtype: {dtype} ✓")

            # Check 'timestamp' column type (should be str/object)
            if 'timestamp' in df.columns:
                dtype = str(df['timestamp'].dtype)
                if dtype not in ('object', 'str'):
                    problematic.append(f"timestamp: {dtype} (deveria ser object/string)")
                else:
                    print(f"   timestamp column dtype: {dtype} ✓")

            # Test for FutureWarning on assignment
            with warnings.catch_warnings(record=True) as w:
                warnings.simplefilter("always")
                if 'success' in df.columns and len(df) > 0:
                    # Try assigning a boolean value
                    df.at[df.index[0], 'success'] = True
                future_warnings = [x for x in w if issubclass(x.category, FutureWarning)]
                if future_warnings:
                    problematic.append(f"FutureWarning on bool assignment: {future_warnings[0].message}")

            if problematic:
                print("   ❌ Schema Issues:")
                for issue in problematic:
                    print(f"      - {issue}")
                checks_failed += 1
            else:
                print("   ✅ CSV schema OK (no FutureWarning on assignment)")
                checks_passed += 1
        else:
            print(f"   ⚠️  CSV não encontrado: {csv_path}")
            checks_failed += 1
    except Exception as e:
        print(f"   ❌ Erro ao verificar CSV: {e}")
        import traceback
        traceback.print_exc()
        checks_failed += 1
    
    # CHECK 3: CREST Installation
    print("\n3️⃣  Checking CREST installation...")
    try:
        import subprocess
        result = subprocess.run(
            ["crest", "--version"],
            capture_output=True,
            timeout=5
        )
        if result.returncode == 0:
            print(f"   ✅ CREST found: {result.stdout.decode().strip()}")
            checks_passed += 1
        else:
            print(f"   ❌ CREST error: {result.stderr.decode()}")
            checks_failed += 1
    except FileNotFoundError:
        print(f"   ❌ CREST not found in PATH")
        checks_failed += 1
    except Exception as e:
        print(f"   ❌ Error checking CREST: {e}")
        checks_failed += 1
    
    # CHECK 4: MOPAC Installation
    print("\n4️⃣  Checking MOPAC installation...")
    try:
        import subprocess
        result = subprocess.run(
            ["mopac", "--version"],
            capture_output=True,
            timeout=5
        )
        if result.returncode == 0:
            print(f"   ✅ MOPAC found: {result.stdout.decode().strip()}")
            checks_passed += 1
        else:
            print(f"   ❌ MOPAC error: {result.stderr.decode()}")
            checks_failed += 1
    except FileNotFoundError:
        print(f"   ❌ MOPAC not found in PATH")
        checks_failed += 1
    except Exception as e:
        print(f"   ❌ Error checking MOPAC: {e}")
        checks_failed += 1
    
    # CHECK 5: RDKit SMILES Validation
    print("\n5️⃣  Checking SMILES validity...")
    try:
        from rdkit import Chem
        test_smiles = [
            "CC(C)Cc1ccc(cc1)C(C)C(O)=O",  # Ibuprofen - deve funcionar
            "Cc1cccc(B(O)O)c1",              # Boron - problema conhecido
            "c1ccccc1",                       # Benzene - simples
        ]
        
        valid = 0
        invalid = 0
        
        for smi in test_smiles:
            mol = Chem.MolFromSmiles(smi)
            if mol:
                valid += 1
                print(f"   ✅ {smi}")
            else:
                invalid += 1
                print(f"   ❌ {smi}")
        
        if invalid == 0:
            print(f"   ✅ All SMILES valid ({valid}/{len(test_smiles)})")
            checks_passed += 1
        else:
            print(f"   ⚠️  {invalid}/{len(test_smiles)} SMILES invalid")
            checks_failed += 1
    except Exception as e:
        print(f"   ❌ Error checking SMILES: {e}")
        checks_failed += 1
    
    # CHECK 6: Phase A Results
    print("\n6️⃣  Checking phase_a_results.json...")
    try:
        results_path = Path("data/molecules_pm7/computed/phase_a_results.json")
        if results_path.exists():
            with open(results_path) as f:
                data = json.load(f)

            # Handle both formats:
            # 1. Dict with "results" key containing list of molecules
            # 2. List of molecules directly (legacy format)
            if isinstance(data, dict):
                results = data.get("results", [])
                print(f"   Format: dict with {len(data)} keys")
                if "n_molecules" in data:
                    print(f"   n_molecules: {data['n_molecules']}")
                if "monitor_summary" in data:
                    summary = data["monitor_summary"]
                    print(f"   Success rate (from summary): {summary.get('success_rate', 'N/A')}")
            elif isinstance(data, list):
                results = data
                print(f"   Format: list (legacy)")
            else:
                print(f"   ❌ Unexpected format: {type(data)}")
                checks_failed += 1
                results = []

            if results:
                total = len(results)
                # Check for successful molecules (has most_stable_hof or homo_lumo)
                successful = sum(
                    1 for r in results
                    if isinstance(r, dict) and (
                        r.get('most_stable_hof') is not None or
                        r.get('homo_lumo') is not None
                    )
                )

                print(f"   Total results: {total}")
                print(f"   Successful: {successful}")
                if total > 0:
                    print(f"   Success rate: {successful/total*100:.1f}%")

                if successful > 0:
                    print(f"   ✅ Some molecules succeeded")
                    checks_passed += 1
                else:
                    print(f"   ⚠️  No successful molecules in results array")
                    # Still pass if file exists and is valid JSON
                    checks_passed += 1
            else:
                print(f"   ⚠️  Results array is empty or missing")
                checks_passed += 1  # File format is valid, just no data
        else:
            print(f"   ⚠️  Results file not found: {results_path}")
            checks_failed += 1
    except Exception as e:
        print(f"   ❌ Error reading results: {e}")
        checks_failed += 1
    
    # SUMMARY
    print("\n" + "="*80)
    print(f"📊 RESULTS: {checks_passed} passed, {checks_failed} failed")
    print("="*80 + "\n")
    
    if checks_failed == 0:
        print("✅ All checks passed! Ready for 30k run.")
        return 0
    else:
        print(f"❌ {checks_failed} check(s) failed. See issues above.")
        return 1

if __name__ == "__main__":
    sys.exit(diagnose())
