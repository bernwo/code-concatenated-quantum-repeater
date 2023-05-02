# Code concatenated Quantum Repeater

Code used to obtain the data for the published paper


Determining the reencoding error probability
--------------------------------------------

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
A version of the data that has been manually reorganized, and also contains columns detailing the ratio between the physical error probability and the reencoding error probability, can be found in `processed_data_for_paper_reorganized.odt`.
The raw data for the paper can be found using the DOI 10.4121/b9c7327e-97b2-4ea2-9b74-18c51f265027.
The data was produced using `tree_channel_error_probability_parameters.csv` as included here as input csv file.
The parameter values for which the simulations were performed as based on the cost-function optimization included in the paper.

For more details about determining the reencoding error probability as a function of the physical error probability, see Supplementary Section S3 of the paper.

