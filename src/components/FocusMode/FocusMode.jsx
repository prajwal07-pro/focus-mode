import React from 'react';
import { useFocusMode } from '../../hooks/useFocusMode';
import FocusModeOverlay from './FocusModeOverlay';
//import styles from './FocusMode.module.css';

const FocusMode = ({ children }) => {
  const { isActive, isScreenSharing, toggleFocusMode, startScreenShare } = useFocusMode();

  return (
    <div className={`${styles.container} ${isActive ? styles.focusActive : ''}`}>
      {children}
      
      {isActive && (
        <FocusModeOverlay 
          onExit={toggleFocusMode}
          isScreenSharing={isScreenSharing}
        />
      )}
      
      <div className={styles.controls}>
        <button onClick={startScreenShare} disabled={isScreenSharing}>
          {isScreenSharing ? 'Screen Sharing Active' : 'Start Screen Share'}
        </button>
        <button onClick={toggleFocusMode}>
          {isActive ? 'Exit Focus Mode' : 'Enter Focus Mode'}
        </button>
      </div>
    </div>
  );
};

export default FocusMode;
