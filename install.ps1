<#
.SYNOPSIS
    rswfn installer for Windows (Windows PowerShell 5.1 and PowerShell 7+).

.DESCRIPTION
    Installs the prebuilt rswfn.exe from GitHub Releases. Needs no Rust
    toolchain and no administrator privileges.

    Typical use:
        irm https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.ps1 | iex

    This script never calls `exit`: the `irm ... | iex` one-liner must not
    close the user's terminal. All work happens inside Install-Rswfn, which
    throws on failure (so `pwsh -File install.ps1` returns a non-zero exit
    code) and returns silently on success.
#>
[CmdletBinding()]
param(
    [Alias('h')]
    [switch]$Help,
    [string]$Version,
    [string]$Dir,
    [string]$BaseUrl,
    [string]$Target,
    [switch]$SkipVerify,
    [switch]$Uninstall,
    [switch]$NoPath
)

# ---------------------------------------------------------------------------
# Output helpers. Progress lines always start with "==>".
# ---------------------------------------------------------------------------
function Write-Step {
    param([string]$Message)
    Write-Host "==> $Message"
}

function Write-WarningLine {
    param([string]$Message)
    Write-Host "==> WARNING: $Message"
}

function Fail {
    param([string]$Message)
    $text = "==> ERROR: $Message"
    Write-Error -Message $text -ErrorAction Continue
    throw $text
}

function Show-Usage {
    Write-Host @'
Usage: install.ps1 [options]

Install the rswfn command-line tool from GitHub Releases.

Options:
  -h, -Help                Show this help and exit.
  -Version <version>       Version to install, e.g. 0.1.0 or v0.1.0.
  -Dir <path>              Install directory (default: %LOCALAPPDATA%\Programs\rswfn).
  -BaseUrl <url>           Download base URL or local release directory.
  -SkipVerify              Do not verify the SHA256 checksum.
  -Target <triple>         Override the detected target (default: x86_64-pc-windows-msvc).
  -Uninstall               Remove rswfn.exe and exit.
  -NoPath                  Do not add the install directory to the user PATH.

Environment variables (command-line options take precedence):
  RSWFN_VERSION, RSWFN_INSTALL_DIR, RSWFN_BASE_URL,
  RSWFN_SKIP_VERIFY, RSWFN_TARGET

Examples:
  irm https://raw.githubusercontent.com/qmlCoder/rswfn/main/install.ps1 | iex
  .\install.ps1 -Version 0.1.0 -NoPath
  .\install.ps1 -Uninstall
'@
}

# ---------------------------------------------------------------------------
# Small helpers.
# ---------------------------------------------------------------------------
function Get-DefaultInstallDir {
    param([switch]$Quiet)
    $localAppData = $env:LOCALAPPDATA
    if ([string]::IsNullOrWhiteSpace($localAppData)) {
        $userProfile = $env:USERPROFILE
        if (-not [string]::IsNullOrWhiteSpace($userProfile)) {
            $localAppData = Join-Path $userProfile 'AppData\Local'
        }
    }
    if ([string]::IsNullOrWhiteSpace($localAppData)) {
        if ($Quiet) { return $null }
        Fail 'LOCALAPPDATA is not set: pass -Dir <path> instead'
    }
    return (Join-Path $localAppData 'Programs\rswfn')
}

function Resolve-InstallDir {
    param(
        [string]$DirPath,
        [switch]$Create
    )
    if ([string]::IsNullOrWhiteSpace($DirPath)) {
        $resolved = Get-DefaultInstallDir
    } else {
        try {
            $resolved = $ExecutionContext.SessionState.Path.GetUnresolvedProviderPathFromPSPath($DirPath)
        } catch {
            Fail "invalid install directory '$DirPath': check the path or pass -Dir <path>"
        }
    }
    if ($Create) {
        if (-not (Test-Path -LiteralPath $resolved -PathType Container)) {
            try {
                New-Item -ItemType Directory -Force -Path $resolved -ErrorAction Stop | Out-Null
            } catch {
                Fail "cannot create install directory '$resolved': check the path or pass -Dir <path>"
            }
        }
    }
    return $resolved
}

function Normalize-PathEntry {
    param([string]$Entry)
    if ($null -eq $Entry) { return '' }
    return $Entry.Trim().TrimEnd('\', '/')
}

function Get-UserPathEntries {
    $userPath = [Environment]::GetEnvironmentVariable('Path', 'User')
    if ([string]::IsNullOrWhiteSpace($userPath)) { return @() }
    return @($userPath.Split(';'))
}

function Add-UserPathEntry {
    param([string]$Directory)
    $entries = Get-UserPathEntries
    $normalizedTarget = Normalize-PathEntry $Directory
    foreach ($entry in $entries) {
        if ((Normalize-PathEntry $entry) -ieq $normalizedTarget) {
            Write-Step "rswfn is already in the user PATH: $Directory"
            return
        }
    }
    $userPath = [Environment]::GetEnvironmentVariable('Path', 'User')
    if ([string]::IsNullOrWhiteSpace($userPath)) {
        $newPath = $Directory
    } else {
        $newPath = $userPath.TrimEnd(';') + ';' + $Directory
    }
    try {
        [Environment]::SetEnvironmentVariable('Path', $newPath, 'User')
    } catch {
        Write-WarningLine "could not update the user PATH automatically: $($_.Exception.Message)"
        Write-WarningLine "add this directory to the PATH manually, then open a new terminal: $Directory"
        return
    }
    Write-Step "added $Directory to the user PATH"
    Write-Step 'restart your terminal (or open a new one) to use the updated PATH'
}

function Remove-UserPathEntry {
    param([string]$Directory)
    $entries = Get-UserPathEntries
    $normalizedTarget = Normalize-PathEntry $Directory
    $kept = @()
    $removed = $false
    foreach ($entry in $entries) {
        if ($entry -eq '') { continue }
        if ((Normalize-PathEntry $entry) -ieq $normalizedTarget) {
            $removed = $true
            continue
        }
        $kept += $entry
    }
    if (-not $removed) { return }
    $newPath = ($kept -join ';')
    try {
        [Environment]::SetEnvironmentVariable('Path', $newPath, 'User')
    } catch {
        Write-WarningLine "could not update the user PATH automatically: $($_.Exception.Message)"
        Write-WarningLine "remove this directory from the PATH manually: $Directory"
        return
    }
    Write-Step "removed $Directory from the user PATH"
}

function Save-RemoteFile {
    param(
        [string]$Url,
        [string]$Destination
    )
    try {
        Invoke-WebRequest -UseBasicParsing -Uri $Url -OutFile $Destination -ErrorAction Stop
        return $null
    } catch {
        return $_
    }
}

function Get-HttpStatus {
    param($ErrorRecord)
    try {
        $response = $ErrorRecord.Exception.Response
        if ($null -ne $response) {
            $statusCode = $response.StatusCode
            if ($null -ne $statusCode) { return [int]$statusCode }
        }
    } catch {
    }
    return 0
}

function Copy-LocalFile {
    param(
        [string]$Source,
        [string]$Destination
    )
    if (-not (Test-Path -LiteralPath $Source -PathType Leaf)) {
        return $null
    }
    try {
        Copy-Item -LiteralPath $Source -Destination $Destination -Force -ErrorAction Stop
        return $true
    } catch {
        return $false
    }
}

function Join-ReleasePath {
    param(
        [string]$Root,
        [string[]]$Parts
    )
    $current = $Root
    foreach ($part in $Parts) {
        $current = Join-Path $current $part
    }
    return $current
}

# ---------------------------------------------------------------------------
# The whole installer lives in this function so that preference variables and
# helper state never leak into the caller's session (irm | iex).
# ---------------------------------------------------------------------------
function Install-Rswfn {
    param(
        [switch]$Help,
        [string]$Version,
        [string]$Dir,
        [string]$BaseUrl,
        [string]$Target,
        [switch]$SkipVerify,
        [switch]$Uninstall,
        [switch]$NoPath
    )

    # These stay local to this function and are gone when it returns or throws.
    $ProgressPreference = 'SilentlyContinue'
    $ErrorActionPreference = 'Continue'

    if ($Help) {
        Show-Usage
        return
    }

    # Environment variables are the default for the matching options.
    if ([string]::IsNullOrWhiteSpace($Version) -and $env:RSWFN_VERSION) { $Version = $env:RSWFN_VERSION }
    if ([string]::IsNullOrWhiteSpace($Dir) -and $env:RSWFN_INSTALL_DIR) { $Dir = $env:RSWFN_INSTALL_DIR }
    if ([string]::IsNullOrWhiteSpace($BaseUrl) -and $env:RSWFN_BASE_URL) { $BaseUrl = $env:RSWFN_BASE_URL }
    if ([string]::IsNullOrWhiteSpace($Target) -and $env:RSWFN_TARGET) { $Target = $env:RSWFN_TARGET }
    if (-not $SkipVerify -and $env:RSWFN_SKIP_VERIFY) {
        switch ($env:RSWFN_SKIP_VERIFY) {
            '1' { $SkipVerify = $true }
            'true' { $SkipVerify = $true }
            'TRUE' { $SkipVerify = $true }
            'yes' { $SkipVerify = $true }
            'YES' { $SkipVerify = $true }
            'on' { $SkipVerify = $true }
            default { }
        }
    }

    # ------------------------------------------------------------------
    # Uninstall: no network, no version/target resolution.
    # ------------------------------------------------------------------
    if ($Uninstall) {
        $installDir = Resolve-InstallDir -DirPath $Dir
        $exePath = Join-Path $installDir 'rswfn.exe'
        if (Test-Path -LiteralPath $exePath -PathType Leaf) {
            try {
                Remove-Item -LiteralPath $exePath -Force -ErrorAction Stop
            } catch {
                Fail "failed to remove '$exePath': close any running rswfn.exe and retry, or pass -Dir <path>"
            }
            Write-Step "removed $exePath"
        } else {
            Write-Step "rswfn is not installed in $installDir; nothing to remove"
        }
        if (Test-Path -LiteralPath $installDir -PathType Container) {
            $remaining = @(Get-ChildItem -LiteralPath $installDir -Force -ErrorAction SilentlyContinue)
            if ($remaining.Count -eq 0) {
                try {
                    Remove-Item -LiteralPath $installDir -Force -ErrorAction Stop
                    Write-Step "removed empty directory $installDir"
                } catch {
                    Write-WarningLine "could not remove empty directory '$installDir': $($_.Exception.Message)"
                }
            }
        }
        $defaultDir = Get-DefaultInstallDir -Quiet
        if (-not [string]::IsNullOrWhiteSpace($defaultDir)) {
            if ((Normalize-PathEntry $installDir) -ieq (Normalize-PathEntry $defaultDir)) {
                Remove-UserPathEntry -Directory $installDir
            }
        }
        return
    }

    # ------------------------------------------------------------------
    # Download base: https:// / http:// / file:/// / plain local directory.
    # ------------------------------------------------------------------
    if ([string]::IsNullOrWhiteSpace($BaseUrl)) {
        $BaseUrl = 'https://github.com/qmlCoder/rswfn/releases/download'
    }
    $baseMode = 'remote'
    $remoteBase = ''
    $localBase = ''
    if ($BaseUrl -match '^https?://') {
        $baseMode = 'remote'
        $remoteBase = $BaseUrl.TrimEnd('/')
        if ([string]::IsNullOrWhiteSpace($remoteBase)) {
            Fail "invalid -BaseUrl '$BaseUrl'"
        }
    } elseif ($BaseUrl -match '^file://') {
        $baseMode = 'local'
        try {
            $localBase = ([Uri]$BaseUrl).LocalPath
        } catch {
            Fail "invalid -BaseUrl '$BaseUrl': expected a file:/// URL or a local directory path"
        }
    } else {
        $baseMode = 'local'
        $localBase = $BaseUrl
    }

    # ------------------------------------------------------------------
    # Version: command line > environment > GitHub releases API.
    # ------------------------------------------------------------------
    if ([string]::IsNullOrWhiteSpace($Version)) {
        if ($baseMode -eq 'local') {
            Fail "no version specified: pass -Version x.y.z when using a local base directory '$localBase' (layout: <base>\v<version>\<asset>)"
        }
        Write-Step 'resolving the latest rswfn release from GitHub'
        try {
            $apiHeaders = @{
                'User-Agent' = 'rswfn-installer'
                'Accept'     = 'application/vnd.github+json'
            }
            $response = Invoke-WebRequest -UseBasicParsing -Uri 'https://api.github.com/repos/qmlCoder/rswfn/releases/latest' -Headers $apiHeaders -ErrorAction Stop
            $releaseInfo = $response.Content | ConvertFrom-Json
            $tagName = $releaseInfo.tag_name
        } catch {
            Fail "could not query the latest release from GitHub: pass -Version x.y.z explicitly (for example: .\install.ps1 -Version 0.1.0). Details: $($_.Exception.Message)"
        }
        if ([string]::IsNullOrWhiteSpace($tagName)) {
            Fail "could not parse 'tag_name' from the GitHub API response: pass -Version x.y.z explicitly"
        }
        $Version = "$tagName"
    }
    $Version = "$Version".Trim()
    if ($Version.StartsWith('v') -or $Version.StartsWith('V')) {
        $Version = $Version.Substring(1)
    }
    if ($Version -notmatch '^[0-9A-Za-z._-]+$') {
        Fail "invalid version '$Version': only characters [0-9A-Za-z._-] are allowed"
    }

    # ------------------------------------------------------------------
    # Target triple.
    # ------------------------------------------------------------------
    if ([string]::IsNullOrWhiteSpace($Target)) {
        $arch = $env:PROCESSOR_ARCHITECTURE
        $nativeArch = $env:PROCESSOR_ARCHITEW6432
        if ($arch -eq 'AMD64') {
            $Target = 'x86_64-pc-windows-msvc'
        } elseif ($arch -eq 'ARM64') {
            Fail 'Windows on ARM64 is not supported yet: no prebuilt rswfn.exe for ARM64 is published; build from source with cargo install rswfn'
        } elseif ($arch -eq 'x86' -and $nativeArch -eq 'AMD64') {
            # 32-bit PowerShell running on 64-bit Windows: the OS is still
            # x86_64, so install the x86_64 binary.
            $Target = 'x86_64-pc-windows-msvc'
        } else {
            Fail "unsupported PROCESSOR_ARCHITECTURE '$arch' (PROCESSOR_ARCHITEW6432='$nativeArch'): pass -Target x86_64-pc-windows-msvc explicitly"
        }
    }
    if ($Target -notmatch '^[0-9A-Za-z._-]+$') {
        Fail "invalid target '$Target': only characters [0-9A-Za-z._-] are allowed"
    }

    # ------------------------------------------------------------------
    # Install directory and temporary working directory.
    # ------------------------------------------------------------------
    $installDir = Resolve-InstallDir -DirPath $Dir -Create
    $targetExe = Join-Path $installDir 'rswfn.exe'

    $tempDir = Join-Path ([IO.Path]::GetTempPath()) ('rswfn-install-' + [Guid]::NewGuid().ToString('N'))
    try {
        New-Item -ItemType Directory -Force -Path $tempDir -ErrorAction Stop | Out-Null
    } catch {
        Fail "cannot create temporary directory '$tempDir': $($_.Exception.Message)"
    }

    try {
        # TLS 1.2 (needed by Windows PowerShell 5.1) and a silent progress bar
        # are set inside this function scope only.
        try {
            [Net.ServicePointManager]::SecurityProtocol = [Net.ServicePointManager]::SecurityProtocol -bor [Net.SecurityProtocolType]::Tls12
        } catch {
        }

        $asset = "rswfn-$Version-$Target.zip"
        $archivePath = Join-Path $tempDir $asset
        $sumsPath = Join-Path $tempDir 'SHA256SUMS'

        # --------------------------------------------------------------
        # Download the archive.
        # --------------------------------------------------------------
        if ($baseMode -eq 'remote') {
            $archiveUrl = "$remoteBase/v$Version/$asset"
            Write-Step "downloading $asset for version $Version"
            $downloadError = Save-RemoteFile -Url $archiveUrl -Destination $archivePath
            if ($null -ne $downloadError) {
                $status = Get-HttpStatus -ErrorRecord $downloadError
                if ($status -eq 404) {
                    Fail "no release asset '$asset' for target '$Target' at $archiveUrl (HTTP 404): check -Version '$Version' or pass -Target x86_64-pc-windows-msvc"
                }
                Fail "failed to download '$archiveUrl': $($downloadError.Exception.Message). Check your network connection or pass -BaseUrl <url>"
            }
        } else {
            $archiveSource = Join-ReleasePath -Root $localBase -Parts @("v$Version", $asset)
            Write-Step "copying $asset from local directory $localBase"
            $copied = Copy-LocalFile -Source $archiveSource -Destination $archivePath
            if ($null -eq $copied -or $copied -eq $false) {
                Fail "archive not found: '$archiveSource'. Check -Version '$Version', -Target '$Target' and the local base directory '$localBase' (layout: <base>\v<version>\<asset>)"
            }
        }

        # --------------------------------------------------------------
        # SHA256 verification.
        # --------------------------------------------------------------
        $haveSums = $false
        if (-not $SkipVerify) {
            if ($baseMode -eq 'remote') {
                $sumsUrl = "$remoteBase/v$Version/SHA256SUMS"
                $sumsError = Save-RemoteFile -Url $sumsUrl -Destination $sumsPath
                $haveSums = ($null -eq $sumsError)
            } else {
                $sumsSource = Join-ReleasePath -Root $localBase -Parts @("v$Version", 'SHA256SUMS')
                $copiedSums = Copy-LocalFile -Source $sumsSource -Destination $sumsPath
                $haveSums = ($copiedSums -eq $true)
            }
            if (-not $haveSums) {
                Write-WarningLine "SHA256SUMS is not available for v${Version}: skipping checksum verification (older release?)"
            }
        } else {
            Write-Step 'skipping checksum verification (-SkipVerify)'
        }

        if ($haveSums -and -not $SkipVerify) {
            $expectedHash = $null
            try {
                $sumsLines = Get-Content -LiteralPath $sumsPath -ErrorAction Stop
            } catch {
                Fail "failed to read '$sumsPath': $($_.Exception.Message)"
            }
            foreach ($line in $sumsLines) {
                if ($line -match '^([0-9A-Fa-f]{64})\s+\*?(.+?)\s*$') {
                    if ($Matches[2] -eq $asset) {
                        $expectedHash = $Matches[1]
                        break
                    }
                }
            }
            if ([string]::IsNullOrWhiteSpace($expectedHash)) {
                Fail "SHA256SUMS does not contain an entry for '$asset': refusing to install an unverified archive (use -SkipVerify only if you trust the download)"
            }
            try {
                $actualHash = (Get-FileHash -LiteralPath $archivePath -Algorithm SHA256 -ErrorAction Stop).Hash
            } catch {
                Fail "failed to compute the SHA256 checksum of '$archivePath': $($_.Exception.Message)"
            }
            if (-not ($expectedHash -ieq $actualHash)) {
                Write-Host "==> ERROR: checksum mismatch for '$asset'"
                Write-Host "==> ERROR:   expected: $($expectedHash.ToLowerInvariant())"
                Write-Host "==> ERROR:   actual:   $($actualHash.ToLowerInvariant())"
                throw "==> ERROR: checksum mismatch for '$asset': nothing was installed. The download may be corrupted or tampered with; retry the download or pass -SkipVerify to bypass this check."
            }
            Write-Step "verified SHA256 checksum of $asset"
        }

        # --------------------------------------------------------------
        # Extract and locate rswfn.exe.
        # --------------------------------------------------------------
        $extractDir = Join-Path $tempDir 'extract'
        try {
            New-Item -ItemType Directory -Force -Path $extractDir -ErrorAction Stop | Out-Null
            Expand-Archive -LiteralPath $archivePath -DestinationPath $extractDir -Force -ErrorAction Stop
        } catch {
            Fail "failed to extract '$asset': $($_.Exception.Message)"
        }

        $candidateExe = Join-Path $extractDir 'rswfn.exe'
        if (-not (Test-Path -LiteralPath $candidateExe -PathType Leaf)) {
            $subDirs = @(Get-ChildItem -LiteralPath $extractDir -Directory -ErrorAction SilentlyContinue)
            if ($subDirs.Count -eq 1) {
                $candidateExe = Join-Path $subDirs[0].FullName 'rswfn.exe'
            }
        }
        if (-not (Test-Path -LiteralPath $candidateExe -PathType Leaf)) {
            Fail "could not find 'rswfn.exe' in '$asset' (unexpected archive layout); please report this at https://github.com/qmlCoder/rswfn/issues"
        }

        # --------------------------------------------------------------
        # Install: replace the existing binary in place.
        # --------------------------------------------------------------
        $stagedExe = Join-Path $installDir ('.rswfn.new.' + [Guid]::NewGuid().ToString('N') + '.exe')
        try {
            Copy-Item -LiteralPath $candidateExe -Destination $stagedExe -Force -ErrorAction Stop
        } catch {
            Fail "failed to copy the binary into '$installDir': check permissions or pass -Dir <path>"
        }
        try {
            Move-Item -LiteralPath $stagedExe -Destination $targetExe -Force -ErrorAction Stop
        } catch {
            Remove-Item -LiteralPath $stagedExe -Force -ErrorAction SilentlyContinue
            Fail "failed to install '$targetExe': close any running rswfn.exe and retry, or pass -Dir <path>"
        }

        # --------------------------------------------------------------
        # Post-install self test: a failure only warns. Older rswfn
        # binaries do not understand --version, and that must not fail
        # the installation.
        # --------------------------------------------------------------
        $versionOutput = ''
        $versionOk = $false
        try {
            $versionOutput = (& $targetExe --version 2>&1 | Out-String)
            $versionOk = ($LASTEXITCODE -eq 0) -and (-not [string]::IsNullOrWhiteSpace($versionOutput))
        } catch {
            $versionOk = $false
        }
        if ($versionOk) {
            Write-Step "self-test ok: $($versionOutput.Trim())"
        } else {
            Write-WarningLine "the installed binary did not report a version for '--version' (normal for very old releases)"
            try {
                $helpOutput = (& $targetExe --help 2>&1 | Out-String)
                if (-not [string]::IsNullOrWhiteSpace($helpOutput)) {
                    Write-Host '==> installed binary help output:'
                    Write-Host $helpOutput.TrimEnd()
                }
            } catch {
            }
        }

        # --------------------------------------------------------------
        # User PATH.
        # --------------------------------------------------------------
        if ($NoPath) {
            Write-Step 'skipping the user PATH update (-NoPath)'
        } else {
            Add-UserPathEntry -Directory $installDir
        }

        Write-Step "rswfn $Version installed to $targetExe"
        Write-Step "run 'rswfn --help' to get started"
    } finally {
        if ($tempDir -and (Test-Path -LiteralPath $tempDir)) {
            try {
                Remove-Item -LiteralPath $tempDir -Recurse -Force -ErrorAction SilentlyContinue
            } catch {
            }
        }
    }
}

Install-Rswfn -Help:$Help -Version "$Version" -Dir "$Dir" -BaseUrl "$BaseUrl" -Target "$Target" -SkipVerify:$SkipVerify -Uninstall:$Uninstall -NoPath:$NoPath
