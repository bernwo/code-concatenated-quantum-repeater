# Code concatenated Quantum Repeater

Code used to obtain the data for the published paper

## 🧮 Calculating effective error probability

_Note: This section applies to the folder `./effective_error_probability/` in this repository._

### Building the effective error probabilities for part of the network

To calculate the secret key rate, which depends on the secret key fraction $f$. This is dependent on the effective error probability of the network, which is what we need to calculate.

This is done via `main.exe`, which is available in `./effective_error_probability/build/bin/` upon building from source from the folder `./effective_error_probability/`. To build from source, run the command (tested on Windows with `MSVC`)

```sh
cmake -S . -B ./build -DCMAKE_BUILD_TYPE=Release; cmake --build ./build --config Release;
```

To generate effective errors for the fault-tolerant error correction protocol with flag qubit, from the `./effective_error_probability/` folder, run

```sh
./build/bin/main.exe flag
```

To generate effective errors for the 1-erasure error correction protocol, from the `./effective_error_probability/` folder, run

```sh
./build/bin/main.exe 1erasure
```

After running `main.exe`, you should see new files being written in `./binary_data/`, which is required by another binary to extrapolate the effective error for the entire network using the approximation derived via a recurrence relation as detailed in the supplemental material. The pre-generated files from running `./secret_key_rate/build/bin/main.exe` are available within this repository in `./effective_error_probability/binary_data/` for your convenience.

## 🔑 Calculating the secret key rate with respect to minimized cost function

Using the data generated in the previous [section](#-calculating-effective-error-probability), we can minimize the cost function and thus calculate the corresponding secret key rate via `./secret_key_rate/build/Release/main.exe`. To build from source, go to `./secret_key_rate/` and run

```sh
cmake -S . -B ./build -DCMAKE_BUILD_TYPE=Release; cmake --build ./build --config Release;
```

For more information, run `./build/Release/main.exe --help`. Running `./build/Release/main.exe` results in data generated in `./binary_data`.

For your convenience, pre-generated data are available in `./secret_key_rate/binary_data` in this repository. To convert these binary data to readable format and easy plotting, we use the Wolfram Mathematica notebook file `./secret_key_rate/post_processing.nb`, which outputs them into `./secret_key_rate/tex/`. Again for your convenience, the pre-converted files are available in `./secret_key_rate/tex/` in this repository.

## 🔍 Determining the re-encoding error probability

_Note: This section applies to the folder `./reencoding_error_probability/` in this repository._

To determine the probability that an error occurs on a logical tree-encoded qubit given that there are loss channels and depolarizing channels on each of the lower-level qubits making it up,
the python file `determine_reencoding_error_probability.py` can be used.
This file depends on the netsquid snippet (a public library built on top of the quantum-network simulator NetSquid) called [NetSquid-TreeCode](https://gitlab.com/softwarequtech/netsquid-snippets/netsquid-treecode).
For instructions on how to use and install it, see the [NetSquid website](https://netsquid.org/snippets/).
The version of the snippet that was used to obtain the data used for the paper is 1.0.1.
Note that while all netsquid snippets are publicly available and open source, installing them (from the netsquid PyPi server) requires netsquid credentials (instructions on how to obtain these can be found on the page of the netsquid website linked above).
These are the same credentials as those that are required to download and install NetSquid itself, or to access the NetSquid forum.

The recommended way to install the snippet is by issuing the command `make requirements` in the terminal when in the folder that holds this README file.
This automatically ensures that the version of the snippet and of NetSquid itself are the same as used to gather the data for the paper.
However, this requires that the user first sets their netsquid credentials as environment variables.
This can be done, for example, by typing the following in the terminal:
`export NETSQUIDPYPI_USER=<netsquid-username>`
`export NETSQUIDPYPI_PWD=<netsquid-password>`
Note that this sets the environment variables only for the current session.
To set them in a persistent manner, the lines can be added to your `~/.bashrc` file (this requires a restart of the terminal to take effect).

When running the python file `determine_reencoding_error_probability.py`, the effective channel on the logical tree-encoded qubit that is induced by loss and depolarizing channels on the lower-level qubits is determined by estimating the Choi state of the channel.
This is done by preparing the state $\tfrac 1 {\sqrt 2} (\ket{00} + \ket{11})$, encoding the second qubit in a tree, subjecting it to loss and noise, and then reencoding it again.
The simulations are performed using NetSquid's stabilizer formalism (see the [NetSquid documentation](https://docs.netsquid.org/latest-release/api_qubits/netsquid.qubits.qformalism.html) for more information).
This formalism allows for efficient simulation of a large number of qubits.
However, the formalism is not able to track mixed states (such as, e.g., the density-matrix formalism).
Therefore, in order obtain a reliable estimate of the Choi state, the simulation needs to be executed many times (we performed every simulation one million times).
By default, multiprocessing is used to speed up the simulation.

What values of the loss, noise and branching vector are used in the simulation is controlled using the `tree_channel_error_probability_parameter.csv` file.
Every row specifies a single combination of parameter values for which the Choi state needs to be determined.
The output of the script is a directory holding all the raw data (one csv file for each value of the parameters), a file holding the combined raw data (`combined_data.csv`), and a file holding the results of performing processing on this data (`processed_data.csv`).
The latter one contains the most interesting and readily-interpreted information.
In this file, each row again corresponds to combination of parameter values, but apart from just listing the parameters it also lists simulation results.
These results include the success probability (i.e., the probability that a tree can actually be decoded after subjected to error and loss), the reencoding error probability (obtained from comparing the fidelity of the Choi state to the fidelity of a Choi state corresponding to a depolarizing channel and determining the depolarizing parameter that results in equal fidelities) and the fidelity to a depolarizing channel (the fidelity between the channel and the depolarizing channel, this is an indicator for how well the channel is approximated by a depolarizing channel).
The processed data also includes standard errors of the estimates and the computation time that was required to obtain the results (in seconds, adding time required by different processes together).
Note that the input csv file does not list the loss probability to which individual qubits are subjected, but instead the parameter `L_0`.
This is the distance over which photonic qubits, assumed to make up the tree code, are transmitted through fiber.
A loss probability is calculated from this assuming an attenuation length of 20 km.
This loss probability is also included in the processed data.

The data that was produced for the paper can be found in `processed_data_for_paper.csv`.
A version of the data that has been manually reorganized, and also contains columns detailing the ratio between the physical error probability and the reencoding error probability, can be found in `processed_data_for_paper_reorganized.ods`.
The raw data for the paper can be found using the DOI 10.4121/b9c7327e-97b2-4ea2-9b74-18c51f265027.
The data was produced using `tree_channel_error_probability_parameters.csv` as included here as input csv file.
The parameter values for which the simulations were performed as based on the cost-function optimization included in the paper.

For more details about determining the reencoding error probability as a function of the physical error probability, see Supplementary Section S3 of the paper.
