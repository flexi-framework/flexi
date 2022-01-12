\hypertarget{installationguidelines}{}

# Installation guidelines \label{chap:installationguidelines}

This chapter contains guidelines to install the code from Github on specific systems.

## Cloning and compiling on machines at the HLRS \label{sec:cloninghlrs}

Unfortunately, the GitHub server is not available on machines at the HLRS, such as the Hazelhen, due to restricted internet access. The workaround is to use ssh tunnels to access the GitHub repositories. Note that the reomte repositories hosted at teh GitLab at the Institute of Aerodynamics and Gasdynamics (IAG), no ssh tunnel is required and cloning works straight forwardly.

The following instructions to access the GitHub repositories on HLRS machines is taken from the HLRS wickie page, see [https://wickie.hlrs.de/platforms/index.php/Secure_Shell_ssh#Git](https://wickie.hlrs.de/platforms/index.php/Secure_Shell_ssh#Git).

### HTTPS

Unfortunately, just using a SSH tunnel as with the SSH and git protocols is not sufficient in this case. Instead, one has to connect via an additional SOCKS proxy on a machine that has unlimited access to the internet, e.g. your local machine.

In order to do so, establish a proxy by using a special feature of OpenSSH:

       ssh -N -D 1080 localhost

This will establish some kind of a "loopback" SSH connection from your local machine to itself which will not execute any command (-N) but act as an SOCKS proxy on port 1080 (-D 1080).

On a second shell, now login to the desired HWW-system and forward a port on the remote machine (e.g. 7777) to the port on your local machine where the newly established SOCKS proxy is listening on (1080):

       ssh -R 7777:localhost:1080 <system-name>.hww.de

By doing so, you have a SOCKS proxy listening on port 7777 of the HWW-system. Hence you can use this proxy for accessing remote git repositories. Unfortunately, the default versions of git installed on the HWW-systems are not capable of doing this. You hence have to load an appropriate version first:

       module load tools/git

In order to use the proxy, you can now add "-c https.proxy='socks5://localhost:7777'" to your git commands, e.g.:

       git -c https.proxy='socks5://localhost:7777' clone https://github.com/flexi-framework/flexi.git

In order to avoid typing this in every git call, you can also set the respective port to be used whenever git talks to a remote repository via HTTPS by

       git config --global https.proxy 'socks5://localhost:7777'

Unfortunately, to connect with GitHub for pulling or pushing, the connection to Hazelhen has to be done via the ssh tunnel.
