#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import logging
import os
import resource

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


def get_memory_bytes(memory_str: str):
    memory_int = int(memory_str.split(" ")[0])
    if memory_str.endswith("GB"):
        return memory_int * 1024 * 1024 * 1024
    elif memory_str.endswith("MB"):
        return memory_int * 1024 * 1024
    elif memory_str.endswith("KB"):
        return memory_int * 1024
    else:
        raise ValueError(f"Invalid memory string: {memory_str}")


def get_soft_limit(value: int):
    return max(1, int(value * 0.9))


def set_max_resources(
    max_cpus: int,
    max_memory_str: str,
    limit_polars: bool = False,
    multiprocess: bool = False,
):
    logger.info(f"Setting max resources to {max_cpus} CPUs and {max_memory_str}")

    max_memory = get_memory_bytes(max_memory_str)
    memory_soft_limit = get_soft_limit(max_memory)
    try:
        # RLIMIT_DATA limits data segment (more reliable than RLIMIT_AS on some systems)
        resource.setrlimit(resource.RLIMIT_DATA, (memory_soft_limit, max_memory))
    except (ValueError, OSError):
        # Fallback to RLIMIT_AS if RLIMIT_DATA doesn't work
        resource.setrlimit(resource.RLIMIT_AS, (memory_soft_limit, max_memory))

    # if running polars, resource.setrlimit returns an error with the cpus
    # instead, we limit the number of threads used by polars
    cpu_soft_limit = get_soft_limit(max_cpus)
    if limit_polars:
        os.environ["POLARS_MAX_THREADS"] = str(cpu_soft_limit)
    elif multiprocess:
        os.environ["OMP_NUM_THREADS"] = str(cpu_soft_limit)
    else:
        resource.setrlimit(resource.RLIMIT_NPROC, (cpu_soft_limit, max_cpus))
