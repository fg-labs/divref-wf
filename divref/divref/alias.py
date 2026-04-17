from typing import Any
from typing import TypeAlias

HailPath: TypeAlias = str
"""Type alias for filesystem paths accepted by Hail: local, GCS (gs://), or HDFS (hdfs://)."""

HailExpression: TypeAlias = Any
"""Type alias for Hail DSL expression objects, which are opaque to the mypy type checker."""
