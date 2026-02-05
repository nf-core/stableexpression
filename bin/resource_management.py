#!/usr/bin/env python3

# Written by Olivier Coen. Released under the MIT license.

import logging
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


def set_max_memory(max_memory_str: str):
    logger.info(f"Setting max memory to {max_memory_str}")

    max_memory = get_memory_bytes(max_memory_str)
    memory_soft_limit = get_soft_limit(max_memory)
    try:
        # RLIMIT_DATA limits data segment (more reliable than RLIMIT_AS on some systems)
        resource.setrlimit(resource.RLIMIT_DATA, (memory_soft_limit, max_memory))
    except (ValueError, OSError):
        # Fallback to RLIMIT_AS if RLIMIT_DATA doesn't work
        resource.setrlimit(resource.RLIMIT_AS, (memory_soft_limit, max_memory))
