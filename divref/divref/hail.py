import os
from pathlib import Path

import hail as hl
import pyspark


def hail_init(
    gcs_credentials_path: Path, spark_driver_memory_gb: int = 1, spark_executor_memory_gb: int = 1
) -> None:
    """
    Initialize Hail with GCS connector and credential configuration.

    Sets ``GOOGLE_APPLICATION_CREDENTIALS`` so the JVM subprocess inherits it,
    then starts Hail with the GCS connector JAR on the Spark classpath.

    Args:
        gcs_credentials_path: Absolute path to a GCP Application Default Credentials
            JSON file. If the file exists and ``GOOGLE_APPLICATION_CREDENTIALS`` is not
            already set, it is exported to the environment before Hail starts.
        spark_driver_memory_gb: Memory in GB to allocate to the Spark driver.
        spark_executor_memory_gb: Memory in GB to allocate to the Spark executor.
    """
    if spark_driver_memory_gb < 1:
        raise ValueError(
            f"Spark driver memory must be at least 1GB. Saw {spark_driver_memory_gb}GB."
        )
    if spark_executor_memory_gb < 1:
        raise ValueError(
            f"Spark driver memory must be at least 1GB. Saw {spark_driver_memory_gb}GB."
        )

    os.environ["PYSPARK_SUBMIT_ARGS"] = (
        f"--driver-memory {spark_driver_memory_gb}g "
        f"--executor-memory {spark_executor_memory_gb}g "
        "pyspark-shell"
    )

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
        }
    )
