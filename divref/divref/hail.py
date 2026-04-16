import os
from pathlib import Path

import hail as hl
import pyspark


def hail_init(gcs_credentials_path: Path) -> None:
    """Initialize Hail with correct configuration."""
    # Set GOOGLE_APPLICATION_CREDENTIALS before hl.init() so the JVM subprocess
    # inherits it. The GCS connector reads this env var to locate ADC credentials;
    # Hadoop-level auth config alone is insufficient when running outside GCP.
    if gcs_credentials_path.exists() and "GOOGLE_APPLICATION_CREDENTIALS" not in os.environ:
        os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = str(gcs_credentials_path)

    jars_dir = os.path.join(pyspark.__path__[0], "jars")
    gcs_jar = os.path.join(jars_dir, "gcs-connector.jar")

    hl.init(
        spark_conf={
            "spark.jars": gcs_jar,
            "spark.driver.extraClassPath": gcs_jar,
            "spark.hadoop.fs.gs.impl": "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFileSystem",
            "spark.hadoop.fs.AbstractFileSystem.gs.impl": "com.google.cloud.hadoop.fs.gcs.GoogleHadoopFS",  # noqa: E501
        }
    )
