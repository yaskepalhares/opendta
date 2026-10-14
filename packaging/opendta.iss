; Inno Setup script for the OpenDTA Windows installer (64-bit, Intel/AMD).
;
;   pyinstaller packaging/opendta.spec --noconfirm
;   iscc /DMyAppVersion=0.2.0 packaging/opendta.iss
;
; Builds dist/OpenDTA-windows-x64-setup.exe from the dist/OpenDTA folder made by
; PyInstaller. tools/build_windows.ps1 runs the whole process. The installer
; works without administrator rights (per-user install); with rights, it can
; install for all users.

#ifndef MyAppVersion
  #define MyAppVersion "0.0.0"
#endif

[Setup]
AppId={{918B8828-E0AB-4114-AC0D-38684306A428}
AppName=OpenDTA
AppVersion={#MyAppVersion}
AppVerName=OpenDTA {#MyAppVersion}
AppPublisher=OpenDTA contributors
AppPublisherURL=https://github.com/yaskepalhares/opendta
AppSupportURL=https://github.com/yaskepalhares/opendta/issues
AppUpdatesURL=https://github.com/yaskepalhares/opendta/releases
DefaultDirName={autopf}\OpenDTA
DefaultGroupName=OpenDTA
DisableProgramGroupPage=yes
PrivilegesRequired=lowest
PrivilegesRequiredOverridesAllowed=dialog
ArchitecturesAllowed=x64compatible
ArchitecturesInstallIn64BitMode=x64compatible
MinVersion=10.0
OutputDir=..\dist
OutputBaseFilename=OpenDTA-windows-x64-setup
SetupIconFile=..\assets\icons\opendta\windows\opendta.ico
UninstallDisplayIcon={app}\OpenDTA.exe
UninstallDisplayName=OpenDTA {#MyAppVersion}
LicenseFile=..\LICENSE
Compression=lzma2/max
SolidCompression=yes
WizardStyle=modern
ChangesAssociations=yes

[Languages]
Name: "english"; MessagesFile: "compiler:Default.isl"
Name: "brazilianportuguese"; MessagesFile: "compiler:Languages\BrazilianPortuguese.isl"

[Tasks]
Name: "desktopicon"; Description: "{cm:CreateDesktopIcon}"; GroupDescription: "{cm:AdditionalIcons}"; Flags: unchecked

[Files]
Source: "..\dist\OpenDTA\*"; DestDir: "{app}"; Flags: ignoreversion recursesubdirs createallsubdirs

[Icons]
Name: "{autoprograms}\OpenDTA"; Filename: "{app}\OpenDTA.exe"
Name: "{autodesktop}\OpenDTA"; Filename: "{app}\OpenDTA.exe"; Tasks: desktopicon

[Registry]
; "Open with" for do-files and .dta datasets; OpenDTA does not become the
; default program for these extensions (it only joins the "Open with" list)
Root: HKA; Subkey: "Software\Classes\OpenDTA.File"; ValueType: string; ValueName: ""; ValueData: "OpenDTA file"; Flags: uninsdeletekey
Root: HKA; Subkey: "Software\Classes\OpenDTA.File\DefaultIcon"; ValueType: string; ValueName: ""; ValueData: "{app}\OpenDTA.exe,0"
Root: HKA; Subkey: "Software\Classes\OpenDTA.File\shell\open\command"; ValueType: string; ValueName: ""; ValueData: """{app}\OpenDTA.exe"" ""%1"""
Root: HKA; Subkey: "Software\Classes\.do\OpenWithProgids"; ValueType: none; ValueName: "OpenDTA.File"; Flags: uninsdeletevalue
Root: HKA; Subkey: "Software\Classes\.ado\OpenWithProgids"; ValueType: none; ValueName: "OpenDTA.File"; Flags: uninsdeletevalue
Root: HKA; Subkey: "Software\Classes\.dta\OpenWithProgids"; ValueType: none; ValueName: "OpenDTA.File"; Flags: uninsdeletevalue
Root: HKA; Subkey: "Software\Classes\Applications\OpenDTA.exe"; ValueType: string; ValueName: "FriendlyAppName"; ValueData: "OpenDTA"; Flags: uninsdeletekey
Root: HKA; Subkey: "Software\Classes\Applications\OpenDTA.exe\SupportedTypes"; ValueType: string; ValueName: ".do"; ValueData: ""
Root: HKA; Subkey: "Software\Classes\Applications\OpenDTA.exe\SupportedTypes"; ValueType: string; ValueName: ".ado"; ValueData: ""
Root: HKA; Subkey: "Software\Classes\Applications\OpenDTA.exe\SupportedTypes"; ValueType: string; ValueName: ".dta"; ValueData: ""

[Run]
Filename: "{app}\OpenDTA.exe"; Description: "{cm:LaunchProgram,OpenDTA}"; Flags: nowait postinstall skipifsilent
