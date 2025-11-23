# Data Directory

This directory contains or will contain observational data for Vela SNR analysis.

## Data Sources

### Primary: Chandra X-ray Observatory
- ObsIDs: 10135, 728, 2202
- Access: https://cxc.cfa.harvard.edu/cda/

### Fallback: Synthetic Clumpy ISM
- Model: White & Long (1991)
- Seed: 42 (reproducible)
- Parameters: 20 random clouds, NE/SW asymmetry 1.28

## File Formats
- FITS files: `.fits` or `.fits.gz`
- Images: `.jpg`, `.png`

## Notes
- All large data files are `.gitignore`d
- Download happens automatically when running script
