from docutils import nodes
from docutils.parsers.rst import Directive


class CommentDirective(Directive):
    has_content = True

    def run(self):
        # Just ignore the content
        return [nodes.comment()]


def setup(app):
    app.add_directive('comment', CommentDirective)
    return {'version': '0.1'}   # returns arbitrary data to signal load success
