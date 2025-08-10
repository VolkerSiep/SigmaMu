# TODO: remove file, but it is used in some tests?
#    use store.predefined_parameters instead.

# from yaml import safe_load
#
# from simu.core.thermo.parameters import StringDictThermoSource
# from simu.app.data import DATA_DIR
#
#
# class PredefinedThermoSource(StringDictThermoSource):
#     def __init__(self):
#         with open(DATA_DIR / "parameters.json") as file:
#             data = safe_load(file)
#         StringDictThermoSource.__init__(self, data)
