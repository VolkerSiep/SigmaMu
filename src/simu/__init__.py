# -*- coding: utf-8 -*-
from logging import getLogger, NullHandler

from ._version import VERSION as __version__
from simu.core.model import Model

logger = getLogger(__name__)
logger.addHandler(NullHandler())
