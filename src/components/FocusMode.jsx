import React, { useEffect, useRef, useState } from 'react';

export default function FocusMode({ sessionId }) {
  const iframeRef = useRef(null);
  const [isScreenSharing, setIsScreenSharing] = useState(false);
  const [screenStream, setScreenStream] = useState(null);
  const [error, setError] = useState('');
  const [status, setStatus] = useState('Ready to share screen');

  useEffect(() => {
    // Listen for predictions from the iframe HTML detector when webcam is active
    const handleMessage = async (event) => {
      if (event.origin !== window.location.origin) return;
      
      if (event.data.type === 'GESTURE_PREDICTION') {
        const { label, probability } = event.data;
        console.log('Focus event:', { label, probability });
      }
    };

    window.addEventListener('message', handleMessage);
    
    return () => {
      window.removeEventListener('message', handleMessage);
      // Clean up screen sharing when component unmounts
      if (screenStream) {
        screenStream.getTracks().forEach(track => track.stop());
      }
    };
  }, [sessionId, screenStream]);

  const startScreenShare = async () => {
    setError('');
    setStatus('Requesting screen access...');
    
    try {
      const stream = await navigator.mediaDevices.getDisplayMedia({
        video: true,
        audio: false
      });

      const [videoTrack] = stream.getVideoTracks();
      const settings = videoTrack.getSettings();

      // Ensure user selected entire screen (not just a window or tab)
      if (settings.displaySurface !== 'monitor') {
        videoTrack.stop();
        setError("Please share your ENTIRE SCREEN. Click the button to try again.");
        setStatus('Screen share failed - wrong selection');
        return;
      }

      setScreenStream(stream);
      setIsScreenSharing(true);
      setStatus('Screen sharing active - Focus detection running');

      // Listen for when user ends screen sharing
      videoTrack.addEventListener('ended', () => {
        console.log("User ended screen share");
        stopScreenShare();
      });

    } catch (error) {
      if (error.name === 'NotAllowedError') {
        setError("Screen share request was cancelled.");
        setStatus('Screen share cancelled');
      } else {
        setError("Error starting screen share: " + error.message);
        setStatus('Screen share error');
      }
      console.error("Error starting screen share:", error);
    }
  };

  const stopScreenShare = () => {
    if (screenStream) {
      screenStream.getTracks().forEach(track => track.stop());
      setScreenStream(null);
    }
    setIsScreenSharing(false);
    setStatus('Screen sharing stopped');
    setError('');
  };

  return (
    <div className="focus-mode-container">
      {/* Screen Share Controls */}
      <div className="screen-share-controls">
        <button 
          onClick={isScreenSharing ? stopScreenShare : startScreenShare}
          className={`share-button ${isScreenSharing ? 'active' : ''}`}
          disabled={status.includes('Requesting')}
        >
          {isScreenSharing ? '🛑 Stop Screen Share' : '📺 Share Your Screen'}
        </button>
        
        <div className="status-info">
          <div className="status-text">{status}</div>
          {error && <div className="error-text">{error}</div>}
        </div>
      </div>

      {/* Main Content Area */}
      <div className="main-content">
        {isScreenSharing ? (
          /* Screen Share Display */
          <div className="screen-share-display">
            <video
              ref={(video) => {
                if (video && screenStream) {
                  video.srcObject = screenStream;
                }
              }}
              autoPlay
              muted
              className="screen-video"
            />
            <div className="screen-overlay">
              <div className="overlay-text">🎯 Analyzing your screen for focus patterns...</div>
              <div className="focus-indicators">
                <div className="indicator">👁️ Eye Tracking: Active</div>
                <div className="indicator">🧠 Focus Analysis: Running</div>
                <div className="indicator">📊 Data Collection: Live</div>
              </div>
            </div>
          </div>
        ) : (
          /* Default Focus Mode with Webcam */
          <div className="webcam-mode">
            <iframe
              ref={iframeRef}
              src="/focus-detector.html"
              style={{
                width: '100%',
                height: '480px',
                border: 'none',
                borderRadius: '8px',
                display: 'block',
                background: '#0b0f1a'
              }}
              title="Focus Detection"
              allow="camera"
            />
            <div className="mode-indicator">
              <div className="small">📹 Webcam Mode - Face focus detection active</div>
            </div>
          </div>
        )}
      </div>

      {/* Usage Instructions */}
      <div className="usage-info">
        <div className="info-section">
          <h4>Screen Share Mode:</h4>
          <ul>
            <li>• Monitors your entire screen for focus and distraction patterns</li>
            <li>• Analyzes window switching, idle time, and productivity apps</li>
            <li>• Provides comprehensive study session insights</li>
          </ul>
        </div>
        <div className="info-section">
          <h4>Webcam Mode:</h4>
          <ul>
            <li>• Face-based focus detection using gesture recognition</li>
            <li>• Monitors attention levels through facial expressions</li>
            <li>• Real-time feedback on engagement and distraction</li>
          </ul>
        </div>
      </div>
    </div>
  );
}
