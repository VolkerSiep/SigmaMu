class HierarchyHandler(dict):
    def __init__(self):
        pass

    def __setitem__(self, name, model):
        if name in self:
            raise KeyError(f"Sub-model of name {name} already defined")
        instance = model.instance()
        super().__setitem__(name, instance)
        return instance