from abc import ABC
from .base import Model


class AModel(Model, ABC):
    """This is a wrapper class to enable short notations for common method
    calls, and by that allowing to shorten the syntax when creating models,
    while still keeping the underlying structure.

    The abbreviations are given in the following table:

    .. currentmodule:: simu.core.model

    ============ ==============================================================
    Abbreviation ... resolves to
    ============ ==============================================================
    pa           :attr:`parameters <simu.Model.parameters>`
    pad          :meth:`parameters.define <parameter.ParameterHandler.define>`
    pas          :meth:`parameters.static <parameter.ParameterHandler.static>`
    pr           :attr:`properties <simu.Model.properties>`
    prd          :meth:`properties.declare <property.PropertyHandler.declare>`
    h            :attr:`hierarchy <simu.Model.hierarchy>`
    hd           :meth:`hierarchy.declare <hierarchy.HierarchyHandler.declare>`
    ha           :meth:`hierarchy.add <hierarchy.HierarchyHandler.add>`
    m            :attr:`materials <simu.Model.materials>`
    md           :meth:`materials.define_port <material.MaterialHandler.define_port>`
    mcf          :meth:`materials.create_flow <material.MaterialHandler.create_flow>`
    mcs          :meth:`materials.create_state <material.MaterialHandler.create_state>`
    r            :attr:`residuals <simu.Model.residuals>`
    ra           :meth:`residuals.add <simu.core.utilities.residual.ResidualHandler.add>`
    ============ ==============================================================
    """

    # parameter handler methods
    @property
    def pa(self):
        return self.parameters

    @property
    def pad(self):
        return self.parameters.define

    @property
    def pas(self):
        return self.parameters.static

    # property handler methods
    @property
    def pr(self):
        return self.properties

    @property
    def prd(self):
        return self.properties.declare

    # hierarchy handler methods
    @property
    def h(self):
        return self.hierarchy

    @property
    def hd(self):
        return self.hierarchy.declare

    def ha(self):
        return self.hierarchy.add

    # materials handler methods
    @property
    def m(self):
        return self.materials

    @property
    def md(self):
        return self.materials.define_port

    @property
    def mcf(self):
        return self.materials.create_flow

    @property
    def mcs(self):
        return self.materials.create_state

    # residuals handler methods
    @property
    def r(self):
        return self.residuals

    @property
    def ra(self):
        return self.residuals.add


