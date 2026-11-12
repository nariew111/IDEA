import os
import subprocess
import logging
import sys
from pathlib import Path

def setup_logger(name, log_file, level=logging.DEBUG):
    """Function to setup as many loggers as you want"""
    log_file = Path(log_file)
    log_file.parent.mkdir(parents=True, exist_ok=True)

    # Clear existing handlers
    logger = logging.getLogger(name)
    if logger.hasHandlers():
        logger.handlers.clear()

    logger.setLevel(level)

    # File handler
    fh = logging.FileHandler(log_file, mode='w')
    fh.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logger.addHandler(fh)

    # Console handler
    ch = logging.StreamHandler(sys.stdout)
    ch.setFormatter(logging.Formatter('%(levelname)s - %(message)s'))
    logger.addHandler(ch)

    return logger

def run_cmd(cmd_list, logger, is_bcftools=False, bcftools_path=None):
    """Runs a command, logs output, and handles errors robustly."""
    cmd_str = ' '.join(map(str, cmd_list))
    logger.info(f"Running command: {cmd_str}")

    env = os.environ.copy()
    if is_bcftools and bcftools_path:
        # Add parent of bcftools to LD_LIBRARY_PATH if not already there
        bcftools_dir = Path(bcftools_path).parent
        ld_path = env.get('LD_LIBRARY_PATH', '')
        if str(bcftools_dir) not in ld_path.split(':'):
            env['LD_LIBRARY_PATH'] = f"{bcftools_dir}:{ld_path}"
            logger.debug(f"Updated LD_LIBRARY_PATH: {env['LD_LIBRARY_PATH']}")

    try:
        # Using subprocess.run for simplicity and better error handling
        result = subprocess.run(
            cmd_list, 
            capture_output=True, 
            text=True, 
            env=env,
            check=False  # We check the return code manually
        )

        if result.stdout:
            logger.debug(f"STDOUT:\n{result.stdout.strip()}")
        if result.stderr:
            # Log stderr as debug, unless return code is non-zero
            if result.returncode == 0:
                logger.debug(f"STDERR:\n{result.stderr.strip()}")
            else:
                logger.error(f"STDERR:\n{result.stderr.strip()}")

        # Raise an exception if the command failed
        if result.returncode != 0:
            raise subprocess.CalledProcessError(
                result.returncode, cmd_str, output=result.stdout, stderr=result.stderr
            )
        
        logger.info(f"Command finished successfully: {cmd_str}")
        return result

    except FileNotFoundError:
        logger.error(f"Command not found: {cmd_list[0]}")
        raise
    except subprocess.CalledProcessError as e:
        logger.error(f"Command failed with exit code {e.returncode}: {cmd_str}")
        # No need to log stderr again as it's already logged above
        raise
    except Exception as e:
        logger.error(f"An unexpected error occurred while running command: {cmd_str}", exc_info=True)
        raise