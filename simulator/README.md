# Try Out Development Containers: C++

[![Open in Dev Containers](https://img.shields.io/static/v1?label=Dev%20Containers&message=Open&color=blue&logo=visualstudiocode)](https://vscode.dev/redirect?url=vscode://ms-vscode-remote.remote-containers/cloneInVolume?url=https://github.com/microsoft/vscode-remote-try-cpp)

A **development container** is a running container with a well-defined tool/runtime stack and its prerequisites. You can try out development containers with **[GitHub Codespaces](https://github.com/features/codespaces)** or **[Visual Studio Code Dev Containers](https://aka.ms/vscode-remote/containers)**.

This is a sample project that lets you try out either option in a few easy steps. We have a variety of other [vscode-remote-try-*](https://github.com/search?q=org%3Amicrosoft+vscode-remote-try-&type=Repositories) sample projects, too.

> **Note:** If you already have a Codespace or dev container, you can jump to the [Things to try](#things-to-try) section.

## Setting up the development container

### GitHub Codespaces
Follow these steps to open this sample in a Codespace:
1. Click the **Code** drop-down menu.
2. Click on the **Codespaces** tab.
3. Click **Create codespace on main** .

For more info, check out the [GitHub documentation](https://docs.github.com/en/free-pro-team@latest/github/developing-online-with-codespaces/creating-a-codespace#creating-a-codespace).

### VS Code Dev Containers

If you already have VS Code and Docker installed, you can click the badge above or [here](https://vscode.dev/redirect?url=vscode://ms-vscode-remote.remote-containers/cloneInVolume?url=https://github.com/microsoft/vscode-remote-try-cpp) to get started. Clicking these links will cause VS Code to automatically install the Dev Containers extension if needed, clone the source code into a container volume, and spin up a dev container for use.

Follow these steps to open this sample in a container using the VS Code Dev Containers extension:

1. If this is your first time using a development container, please ensure your system meets the pre-reqs (i.e. have Docker installed) in the [getting started steps](https://aka.ms/vscode-remote/containers/getting-started).

2. To use this repository, you can either open the repository in an isolated Docker volume:

    - Press <kbd>F1</kbd> and select the **Dev Containers: Try a Sample...** command.
    - Choose the "C++" sample, wait for the container to start, and try things out!
        > **Note:** Under the hood, this will use the **Dev Containers: Clone Repository in Container Volume...** command to clone the source code in a Docker volume instead of the local filesystem. [Volumes](https://docs.docker.com/storage/volumes/) are the preferred mechanism for persisting container data.

   Or open a locally cloned copy of the code:

   - Clone this repository to your local filesystem.
   - Press <kbd>F1</kbd> and select the **Dev Containers: Open Folder in Container...** command.
   - Select the cloned copy of this folder, wait for the container to start, and try things out!

## Things to try

Once you have this sample opened, you'll be able to work with it like you would locally.

Some things to try:

1. **Edit:**
   - Open `main.cpp`
   - Try adding some code and check out the language features.
   - Make a spelling mistake and notice it is detected. The [Code Spell Checker](https://marketplace.visualstudio.com/items?itemName=streetsidesoftware.code-spell-checker) extension was automatically installed because it is referenced in `.devcontainer/devcontainer.json`.
   - Also notice that utilities like `Vcpkg` and the [C++](https://marketplace.visualstudio.com/items?itemName=ms-vscode.cpptools) extension are installed. Tools are installed in the `mcr.microsoft.com/devcontainers/cpp` image and Dev Container settings and metadata are automatically picked up from [image labels](https://containers.dev/implementors/reference/#labels).

1. **Terminal:** Press <kbd>ctrl</kbd>+<kbd>shift</kbd>+<kbd>\`</kbd> and type `uname` and other Linux commands from the terminal window.

1. **Build, Run, and Debug:**
   - Open `main.cpp`
   - Add a breakpoint (e.g. on line 7).
   - Press <kbd>F5</kbd> to launch the app in the container.
   - Once the breakpoint is hit, try hovering over variables, examining locals, and more.

1. **Install the GitHub CLI using a Dev Container Feature:**
   - Press <kbd>F1</kbd> and select the **Dev Containers: Configure Container Features...** or **Codespaces: Configure Container Features...** command.
   - Type "github" in the text box at the top.
   - Check check box next to "GitHub CLI" (published by devcontainers) 
   - Click OK
   - Press <kbd>F1</kbd> and select the **Dev Containers: Rebuild Container** or **Codespaces: Rebuild Container** command so the modifications are picked up.

## Contributing

This project welcomes contributions and suggestions. Most contributions require you to agree to a
Contributor License Agreement (CLA) declaring that you have the right to, and actually do, grant us
the rights to use your contribution. For details, visit https://cla.microsoft.com.

When you submit a pull request, a CLA-bot will automatically determine whether you need to provide
a CLA and decorate the PR appropriately (e.g., label, comment). Simply follow the instructions
provided by the bot. You will only need to do this once across all repos using our CLA.

This project has adopted the [Microsoft Open Source Code of Conduct](https://opensource.microsoft.com/codeofconduct/).
For more information see the [Code of Conduct FAQ](https://opensource.microsoft.com/codeofconduct/faq/) or
contact [opencode@microsoft.com](mailto:opencode@microsoft.com) with any additional questions or comments.

## License

Copyright © Microsoft Corporation All rights reserved.<br />
Licensed under the MIT License. See LICENSE in the project root for license information.




- <kbd>CMD + Shift + B</kbd> to build
- <kbd>CMD + Shift + P</kbd> => Tasks: Run Task => Run Executable


## TODO
- Improvement: consider **cling** for a compiler.
- Improvement: test library, e.g. [gtest](https://matgomes.com/integrate-google-test-into-cmake/) or [Catch2](https://github.com/catchorg/Catch2).
- Improvement: package manager, e.g. **conan**.
- Improvement: strong types or unit libraries in the OrbitalMechanics namespace (**kilometers** and **radians** instead of **double**). Possible solutions: type-safe wrappers using **struct**, type aliases with **std::chrono::duration**, unit libraries like **MPark.Units** or **Boost.Units**.
- Improvement: consider [OTL](https://github.com/Jmbryan/OTL) for orbital mechanics.
- Improvement: refactor TLE into Orbit
- Improvement: precision of ecef2eci conversion.

https://programiz.pro/ide/python
https://raw.githubusercontent.com/eribean/Geneci/refs/heads/main/geneci/conversions.py

# Define the input parameters
ecef_point = np.array([3784972.26595145, 896365.783633651, 5037929.21063147], dtype=np.float64)
utc_time = datetime.strptime("04-Jul-2023 14:25:00", "%d-%b-%Y %H:%M:%S")
ecef_velocity = np.array([0.0, 0.0, 0.0], dtype=np.float64)

# Invoke the function with the given inputs
eci_point, eci_velocity = ecef_to_eci(ecef_point, utc_time, ecef_velocity)

print(eci_point)    # [-3410576.75207984  1849283.01581851  5045625.42664478]
print(eci_velocity) # [-134.8406056  -249.53731679    0.31331844]

## Useful stuffs

[Numerical integration in C++ with Boost odeint](http://boccelliengineering.altervista.org/junk/boost_integration/boost_odeint.html)
