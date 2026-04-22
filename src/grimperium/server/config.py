"""Server configuration via pydantic-settings.

All environment variables use the GRIMPERIUM_ prefix.
"""

from pydantic_settings import BaseSettings, SettingsConfigDict


class ServerConfig(BaseSettings):
    """Configuration for the Grimperium distributed processing server."""

    csv_path: str
    api_token: str = ""
    heartbeat_timeout_s: int = 300
    watchdog_interval_s: int = 30
    startup_grace_s: int = 120
    stuck_assigned_threshold_h: float = 2.0
    host: str = "0.0.0.0"
    port: int = 8000
    max_reruns: int = 3

    model_config = SettingsConfigDict(env_prefix="GRIMPERIUM_")
