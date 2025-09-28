export class KeyboardShortcuts {
  constructor(focusManager) {
    this.focusManager = focusManager;
    this.bindEvents();
  }

  bindEvents() {
    document.addEventListener('keydown', this.handleKeyDown.bind(this));
  }

  handleKeyDown(e) {
    // Ctrl/Cmd + Shift + F for focus mode
    if ((e.ctrlKey || e.metaKey) && e.shiftKey && e.key === 'F') {
      e.preventDefault();
      this.focusManager.toggleFocusMode();
    }
    
    // Ctrl/Cmd + Shift + S for screen share
    if ((e.ctrlKey || e.metaKey) && e.shiftKey && e.key === 'S') {
      e.preventDefault();
      if (!this.focusManager.isScreenSharing()) {
        this.focusManager.startScreenShare();
      }
    }
  }

  destroy() {
    document.removeEventListener('keydown', this.handleKeyDown);
  }
}

export default KeyboardShortcuts;
