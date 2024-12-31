from pathlib import Path
from re import split
from docutils import nodes
from docutils.parsers.rst import Directive
from sphinx.directives.code import LiteralInclude

# Path for example files to be included in documentation
EXAMPLE_PATH = Path(__file__).parents[2] / "src" / "simu" / "examples"

class CommentDirective(Directive):
    has_content = True

    def run(self):
        # Just ignore the content
        return [nodes.comment()]


class ExampleInclude(LiteralInclude):
    """
    Custom directive to include example files with a base path.
    """
    def run(self):
        # ugly MONKEY_PATCH ... let's hope that keeps on working
        # I need to prevent the LiteralInclude.run() method to prepend the
        # path of the document on matter what. Instead, I need to include
        # the file relative to the simu.examples module.

        old_relpath, self.env.relfn2path = self.env.relfn2path, self._relfn2path
        self.docname = self.env.docname
        result = super().run()
        self.env.relfn2path = old_relpath
        return result

    def _relfn2path(self, filename, docname=None):
        """Return paths to a file referenced from a document, relative to
        documentation root and absolute."""
        abs_path = EXAMPLE_PATH / filename
        doc_path = Path(docname or self.docname).parent.absolute()
        rel_path = abs_path.relative_to(doc_path, walk_up=True)
        return str(rel_path), str(abs_path)


def setup(app):
    app.add_directive('comment', CommentDirective)
    app.add_directive("exampleinclude", ExampleInclude)
    return {'version': '0.1'}   # returns arbitrary data to signal load success
