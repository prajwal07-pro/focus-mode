class FocusModeManager {
  constructor() {
    this.mediaStream = null;
    this.focusModeActive = false;
    this.screenShareUI = null;
  }

  async startScreenShare() {
    const controller = new CaptureController();
    
    this.mediaStream = await navigator.mediaDevices.getDisplayMedia({
      controller,
      video: true,
      audio: true
    });

    const [track] = this.mediaStream.getVideoTracks();
    const displaySurface = track.getSettings().displaySurface;
    
    if (displaySurface !== "monitor") {
      controller.setFocusBehavior("focus-capturing-application");
    }

    return this.mediaStream;
  }

  toggleFocusMode() {
    this.focusModeActive = !this.focusModeActive;
    
    if (this.focusModeActive) {
      this.hideFocusModeElements();
    } else {
      this.showFocusModeElements();
    }
  }

  hideFocusModeElements() {
    const elementsToHide = [
      '.navbar', '.sidebar', '.header', '.footer', '#screen-share-controls'
    ];
    
    elementsToHide.forEach(selector => {
      const elements = document.querySelectorAll(selector);
      elements.forEach(el => {
        el.style.display = 'none';
        el.setAttribute('data-hidden-by-focus', 'true');
      });
    });
  }

  showFocusModeElements() {
    const hiddenElements = document.querySelectorAll('[data-hidden-by-focus="true"]');
    hiddenElements.forEach(el => {
      el.style.display = '';
      el.removeAttribute('data-hidden-by-focus');
    });
  }

  isScreenSharing() {
    return this.mediaStream && 
           this.mediaStream.getVideoTracks().length > 0 && 
           this.mediaStream.getVideoTracks()[0].readyState === 'live';
  }
}

export default FocusModeManager;
