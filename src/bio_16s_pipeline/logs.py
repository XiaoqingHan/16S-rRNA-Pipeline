import logging
import multiprocessing
from logging.handlers import QueueHandler, QueueListener
import os
import time

log_queue = None
_log_listener = None

def init_log_queue():
    global log_queue
    if log_queue is None:
        ctx = multiprocessing.get_context("spawn")  # spawn 模式安全
        log_queue = ctx.Queue()
    return log_queue

def _configure_worker_process_logger():
    logger = logging.getLogger('16S_Analysis')
    logger.handlers.clear()
    queue_handler = QueueHandler(log_queue)
    queue_handler.setLevel(logging.INFO)
    logger.addHandler(queue_handler)
    logger.setLevel(logging.INFO)
    logger.propagate = False
    return logger

def init_worker_logger(q):
    global log_queue
    log_queue = q
    _configure_worker_process_logger()

def get_logger():
    return logging.getLogger('16S_Analysis')

def start_log_listener(queue=None, log_filename=None):
    global _log_listener
    global log_queue
    if _log_listener is not None:
        return _log_listener

    if queue is not None:
        log_queue = queue
    if log_queue is None:
        init_log_queue()

    log_dir = "logs"
    os.makedirs(log_dir, exist_ok=True)
    if log_filename is None:
        log_filename = os.path.join(log_dir, f"16s_analysis_{time.strftime('%Y%m%d_%H%M%S')}.log")

    formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - [%(processName)s] - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    file_handler = logging.FileHandler(log_filename, encoding='utf-8')
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)

    _log_listener = QueueListener(log_queue, file_handler, stream_handler)
    _log_listener.start()

    _configure_worker_process_logger()
    return _log_listener

def stop_log_listener():
    global _log_listener
    if _log_listener is not None:
        _log_listener.stop()
        _log_listener = None
    logging.basicConfig(level=logging.INFO)
    get_logger().info("Log listener stopped.")
