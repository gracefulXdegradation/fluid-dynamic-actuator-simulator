[![Open in Dev Containers](https://img.shields.io/static/v1?label=Dev%20Containers&message=Open&color=blue)](https://vscode.dev/redirect?url=vscode://ms-vscode-remote.remote-containers/cloneInVolume?url=https://github.com/gracefulXdegradation/flight-dynamics-simulator)

If you already have VS Code and Docker installed, you can click the badge above or [here](https://vscode.dev/redirect?url=vscode://ms-vscode-remote.remote-containers/cloneInVolume?url=https://github.com/gracefulXdegradation/flight-dynamics-simulator) to get started. Clicking these links will cause VS Code to automatically install the Dev Containers extension if needed, clone the source code into a container volume, and spin up a dev container for use.

# Working with multiple dev containers in VS Code

## Prerequisites

- [Dev Containers Extension](https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.remote-containers) installed in VS Code.

## How to use multiple dev containers in a single VS Code window

- Open local VS Code window on cloned repo.
- From the command palette `Dev Containers: Reopen in Container`, pick Python Container. 
- This will open a new VS Code window connected to the selected container. 
- From the command palette `Dev Containers: Switch Container`, pick Node Container.
- This will reload the VS Code window connected to the selected container.

## How to use multiple dev containers in a multiple VS Code window

1. Open local VS Code window on cloned repo.
2. From the command palette `Dev Containers: Reopen in Container`, pick Python Container.
3. This will open a new VS Code window connected to the selected container.
4. Go to File > New Window.
5. In this new window, from the command palette `Dev Containers: Reopen in Container`, pick Node Container.
6. This will open a new VS Code window connected to the selected container.

## Additional Resources

- [Dev Containers Supporting tools and services](https://containers.dev/supporting)
- [Dev Containers Documentation for VS Code](https://code.visualstudio.com/docs/remote/containers)
- [Dev Containers Documentation](https://containers.dev/)