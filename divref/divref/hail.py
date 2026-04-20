import os
from pathlib import Path

import hail as hl
import pyspark

from divref import defaults


def hail_init(gcs_credentials_path: Path) -> None:
    """
    Initialize Hail with GCS connector and credential configuration.

    Sets ``GOOGLE_APPLICATION_CREDENTIALS`` so the JVM subprocess inherits it,
    then starts Hail with the GCS connector JAR on the Spark classpath.

    Args:
        gcs_credentials_path: Absolute path to a GCP Application Default Credentials
            JSON file. If the file exists and ``GOOGLE_APPLICATION_CREDENTIALS`` is not
            already set, it is exported to the environment before Hail starts.
    """
    if gcs_credentials_path.exists() and "GOOGLE_APPLICATION_CREDENTIALS" not in os.environ:
        os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = str(gcs_credentials_path)

    jars_dir = os.path.join(pyspark.__path__[0], "jars")
    gcs_jar = os.path.join(jars_dir, "gcs-connector.jar")

    if not os.path.exists(gcs_jar):
        raise FileNotFoundError(
            f"GCS connector JAR not found at {gcs_jar}. Run 'pixi run setup-gcs' to download it."
        )

    hl.init(
        spark_conf={
            "spark.jars": gcs_jar,
            "spark.driver.extraClassPath": gcs_jar,
            "spark.hadoop.fs.gs.impl": "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem",
            "spark.hadoop.fs.AbstractFileSystem.gs.impl": "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS",  # noqa: E501
        },
        default_reference=defaults.REFERENCE_GENOME,
    )
