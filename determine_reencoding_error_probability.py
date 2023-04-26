import netsquid as ns
from netsquid_treecode.tree_channel_error_probability import collect_choi_states_from_csv


if __name__ == "__main__":
    ns.set_qstate_formalism(ns.QFormalism.STAB)
    collect_choi_states_from_csv(csv_specification_filename="tree_channel_error_probability_parameters.csv",
                                 number_of_states=1E6, path="data", multiprocessing=True)
