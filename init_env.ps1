$env:PATH = "C:\Program Files\Git\usr\bin;" + $env:PATH
$env:PATH += ";$PSScriptRoot\build-system\scripts"
Write-Host "Environment initialized. 'erbb' command is now available."
function erbb { python "$PSScriptRoot\build-system\scripts\erbb" $args }

