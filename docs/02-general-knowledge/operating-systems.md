---
title: Operating systems
slug: /general-knowledge/os
---

## Shells

A **shell** executes commands you enter in a terminal in order to interact with the operating system. There are various shells that differ not in their overall capabilities but mostly in their syntax, default functions and tools they provide out of the box and ergonomics.
Shells can be used interactively or as a scripting language.

* [**Bash**](https://www.gnu.org/software/bash/manual/bash.html): the default shell on most Linux and older MacOS systems.
* [**zsh**](https://www.zsh.org/): the default shell on MacOS. Very similar to bash in terms of syntax, but not fully compatible. A popular framework for plugins and customizations is [oh my zsh](https://ohmyz.sh/).
* [**fish**](https://fishshell.com/): less widespread but modern shell with a focus on user-friendliness and scripting.
* **sh**: not a shell per se, but a (minimalistic) programming language described by the POSIX standard. Typically `/bin/sh` is a symbolic link to an actual implementation of the language (which *can* be bash, as bash can be run in POSIX compliant mode).
* [**PowerShell**](https://learn.microsoft.com/en-us/powershell/): the default shell on Microsoft Windows.

Find out which shell you are running with the command `echo $0`.

A **shebang line** is the first line of a shell script that indicates the shell in which the script is supposed to be run. It starts with `#!`, followed by the command to run the shell executable. This is `#!/bin/bash` for Bash or `#!/usr/bin/zsh` for zsh.
