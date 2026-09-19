"""Deterministic playback controller shared by the plots and map."""

from __future__ import annotations

from PySide6.QtCore import QObject, QTimer, Qt, Signal


class PlaybackController(QObject):
    """Drive a sample index with a Qt timer while keeping stepping testable."""

    indexChanged = Signal(int)
    playingChanged = Signal(bool)
    speedChanged = Signal(float)

    def __init__(self, sample_count: int = 0, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._sample_count = max(0, int(sample_count))
        self._index = 0
        self._speed = 1.0
        self.timer = QTimer(self)
        self.timer.setTimerType(Qt.TimerType.PreciseTimer)
        self.timer.timeout.connect(self._on_timeout)
        self._update_interval()

    @property
    def sample_count(self) -> int:
        return self._sample_count

    @property
    def index(self) -> int:
        return self._index

    @property
    def speed(self) -> float:
        return self._speed

    @property
    def is_playing(self) -> bool:
        return self.timer.isActive()

    def set_sample_count(self, count: int, *, reset: bool = False) -> None:
        self._sample_count = max(0, int(count))
        if reset:
            self.set_index(0)
        elif self._sample_count:
            self.set_index(min(self._index, self._sample_count - 1))
        else:
            self.set_index(0)
        if not self._sample_count:
            self.pause()

    def set_index(self, index: int) -> None:
        maximum = max(0, self._sample_count - 1)
        bounded = min(max(0, int(index)), maximum)
        if bounded != self._index:
            self._index = bounded
            self.indexChanged.emit(self._index)
        elif self._sample_count == 1:
            # A caller setting the only sample still gets a useful event when
            # it is reloading a lap.
            self.indexChanged.emit(self._index)

    def step(self, count: int = 1) -> int:
        """Advance by ``count`` samples and return the resulting index."""

        if self._sample_count <= 1:
            self.pause()
            return self._index
        target = self._index + int(count)
        if target >= self._sample_count:
            target = self._sample_count - 1
            self.set_index(target)
            self.pause()
        elif target < 0:
            self.set_index(0)
        else:
            self.set_index(target)
        return self._index

    # ``advance`` is a readable alias for controller users that do not care
    # whether the timer or a test initiated the movement.
    advance = step

    def set_speed(self, speed: float) -> None:
        bounded = min(16.0, max(0.05, float(speed)))
        if bounded == self._speed:
            return
        self._speed = bounded
        self._update_interval()
        self.speedChanged.emit(self._speed)

    def play(self) -> None:
        if self._sample_count <= 1:
            return
        if self._index >= self._sample_count - 1:
            self.set_index(0)
        self._update_interval()
        if not self.timer.isActive():
            self.timer.start()
            self.playingChanged.emit(True)

    def pause(self) -> None:
        was_active = self.timer.isActive()
        self.timer.stop()
        if was_active:
            self.playingChanged.emit(False)

    def toggle(self) -> None:
        if self.is_playing:
            self.pause()
        else:
            self.play()

    def _on_timeout(self) -> None:
        self.step()

    def _update_interval(self) -> None:
        # The viewer advances by one source sample at 30 Hz.  The speed slider
        # therefore changes the timer cadence without discarding samples.
        interval_ms = max(1, int(round(1000.0 / (30.0 * self._speed))))
        self.timer.setInterval(interval_ms)
