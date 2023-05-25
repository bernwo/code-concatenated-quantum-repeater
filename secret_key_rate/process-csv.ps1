Get-ChildItem "$env:USERPROFILE\Documents\GitHub\one-way-quantum-repeaters\tex\*\*.csv" | ForEach-Object {(Get-Content $_) -Replace '"Indeterminate"', 'nan' | Set-Content $_};
Get-ChildItem "$env:USERPROFILE\Documents\GitHub\one-way-quantum-repeaters\tex\*points*\*\*.csv" | ForEach-Object {(Get-Content $_) -Replace '"Indeterminate"', 'nan' | Set-Content $_};

Get-ChildItem "$env:USERPROFILE\Documents\GitHub\one-way-quantum-repeaters\tex\mII*\*mII_*.csv" | ForEach-Object {(Get-Content $_) -Replace ',0$', ',nan' | Set-Content $_};
Get-ChildItem "$env:USERPROFILE\Documents\GitHub\one-way-quantum-repeaters\tex\*points*\mII*\*mII_*.csv" | ForEach-Object {(Get-Content $_) -Replace ',0$', ',nan' | Set-Content $_};
Get-ChildItem "$env:USERPROFILE\Documents\GitHub\one-way-quantum-repeaters\tex\mII*\*mII_*.csv" | ForEach-Object {(Get-Content $_) -Replace ',0.$', ',nan' | Set-Content $_};
Get-ChildItem "$env:USERPROFILE\Documents\GitHub\one-way-quantum-repeaters\tex\*points*\mII*\*mII_*.csv" | ForEach-Object {(Get-Content $_) -Replace ',0.$', ',nan' | Set-Content $_};

Get-ChildItem "$env:USERPROFILE\Documents\GitHub\one-way-quantum-repeaters\tex\mI*\*mI_*.csv" | ForEach-Object {(Get-Content $_) -Replace ',0$', ',nan' | Set-Content $_};
Get-ChildItem "$env:USERPROFILE\Documents\GitHub\one-way-quantum-repeaters\tex\*points*\mI*\*mI_*.csv" | ForEach-Object {(Get-Content $_) -Replace ',0$', ',nan' | Set-Content $_};
Get-ChildItem "$env:USERPROFILE\Documents\GitHub\one-way-quantum-repeaters\tex\mI*\*mI_*.csv" | ForEach-Object {(Get-Content $_) -Replace ',0.$', ',nan' | Set-Content $_};
Get-ChildItem "$env:USERPROFILE\Documents\GitHub\one-way-quantum-repeaters\tex\*points*\mI*\*mI_*.csv" | ForEach-Object {(Get-Content $_) -Replace ',0.$', ',nan' | Set-Content $_};

Get-ChildItem "$env:USERPROFILE\Documents\GitHub\one-way-quantum-repeaters\tex\ntot*\*ntot_*.csv" | ForEach-Object {(Get-Content $_) -Replace ',0$', ',nan' | Set-Content $_};
Get-ChildItem "$env:USERPROFILE\Documents\GitHub\one-way-quantum-repeaters\tex\*points*\ntot*\*ntot_*.csv" | ForEach-Object {(Get-Content $_) -Replace ',0$', ',nan' | Set-Content $_};

Write-Host "Done.";