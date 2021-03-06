Hydration free energy calculation workflow:
-------------------------------------------
* Perform COSMO calculation for both solute and solvent to generate the .cosmo files.
   - Parameters in this example were optimized for GEPOL SES surfaces with solvent radius of 1.3 Bohr, BP86 functional and SVP basis set (aug-cc-pvtz was used for H2O).
   - The current .cosmo file format follows the settings in https://github.com/qsimulate/qm/tree/cosmo_file. Other format can be customized in [pycosmosac/cosmo](../../../pycosmosac/cosmo).
* Run [hydration_energy.py](./hydration_energy.py).
