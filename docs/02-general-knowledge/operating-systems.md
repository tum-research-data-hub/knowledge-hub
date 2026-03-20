---
title: Operating systems
slug: /general-knowledge/os
---
# Operating systems
An operating system is the collection of software that manages a computer's resources and acts as a layer between the firmware and application (i.e. user-facing) software.

How you as a user interact with the computer greatly depends on which operating system is running.

## Common operating systems

| Name | Notes | Common usecases | Usage at TUM / LRZ |
| ---- | ---- | ---- | ---- |
| Linux | any operating system that uses the Linux kernel may refer to itself as "Linux", a popular one being "GNU/Linux" | servers, embedded / IoT devices, HPC clusters, personal computers | LRZ compute infrastructure |
| Windows |  | servers, personal computers |  |
| MacOS |  | personal computers |  |
| Android | based on a modified version of the Linux kernel | phones, tablets, embedded / IoT devices |  |
| iOS |  | phones, tablets |  |

### Linux distributions
A Linux distribution is an operating system that is built on top of the Linux kernel and includes a collection of components that shapes the user experience. These components might be a package manager, window manager, init system, network configuration software and further tools and utilities like office software.
While there are [thousands of Linux distributions](https://en.wikipedia.org/wiki/List_of_Linux_distributions), only a few are widespread in certain domains.

| Distribution | Notes | Common usecases | Usage at TUM / LRZ |
| --- | --- | --- | --- |
| Debian | slow release cycle, very stable | servers | webservers ? |
| Ubuntu | derived from Debian | desktop computers, servers |  |
| CentOS |  |  |  |
| AlmaLinux |  |  |  |
| SUSE |  |  |  |
| Fedora |  |  |  |
| Arch | rolling release cycle |  | desktop computers |

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
