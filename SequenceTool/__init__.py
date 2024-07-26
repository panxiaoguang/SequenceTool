from .index import index
from .aiart import aiart
from .chatapp import aichat
from .Blog import blogview
from .Blogs import routers
__all__ = [
    'index',
    'aiart',
    'aichat',
    'blogview'
]

__all__.extend(routers)
