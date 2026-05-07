#!/bin/bash
# scripts/test_phase_a.sh

echo "=== Phase A Integration Test ==="
echo ""

echo "Step 1: Run 3-molecule batch"
python -m grimperium.cli.main <<EOF
1
3
20
20
y
EOF

echo ""
echo "Step 2: Verify CSV exists"
if [ -f "data/molecules_pm7/computed/thermo_pm7.csv" ]; then
    echo "✓ CSV exists"
else
    echo "✗ CSV not found"
    exit 1
fi

echo ""
echo "Step 3: Check column count"
COLS=$(head -n 1 data/molecules_pm7/computed/thermo_pm7.csv | tr ',' '\n' | wc -l)
echo "Columns found: $COLS"
if [ "$COLS" -eq 49 ]; then
    echo "✓ Correct column count (49)"
else
    echo "✗ Wrong column count (expected 49, got $COLS)"
    exit 1
fi

echo ""
echo "Step 4: Run pytest"
pytest tests/test_settings_phase_a.py tests/test_csv_schema.py tests/test_molecule_processor.py -v

echo ""
echo "=== Phase A Test Complete ==="
