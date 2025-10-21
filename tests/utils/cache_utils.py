import logging
from pathlib import Path
from typing import List
from tests.utils.run_ppanggolin import run_ppanggolin_command

logger = logging.getLogger(__name__)


def run_with_cache(request, cache_key: str, outdir: str, cmds: List[str]) -> Path:
    """
    Run a ppanggolin command if no cached result exists, otherwise reuse cache.

    """
    cache = request.config.cache
    cached_dir = cache.get(cache_key, None)
    logger.warning(f"cached_dir {cached_dir}")

    if cached_dir and Path(cached_dir).exists():
        outdir = Path(cached_dir)
        logger.warning(f"Reusing cached result for {cache_key} at {outdir}")

    else:

        for cmd in cmds:
            run_ppanggolin_command(cmd)

        cache.set(cache_key, str(outdir))
        logger.info(f"Cached result for {cache_key} at {outdir}")

    return outdir
