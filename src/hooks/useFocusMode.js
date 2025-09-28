import { useState, useEffect, useRef } from 'react';
import FocusModeManager from '../services/FocusModeManager';

export const useFocusMode = () => {
  const [isActive, setIsActive] = useState(false);
  const [isScreenSharing, setIsScreenSharing] = useState(false);
  const managerRef = useRef(null);

  useEffect(() => {
    managerRef.current = new FocusModeManager();
    
    return () => {
      if (managerRef.current?.mediaStream) {
        managerRef.current.mediaStream.getTracks().forEach(track => track.stop());
      }
    };
  }, []);

  const toggleFocusMode = () => {
    if (managerRef.current) {
      managerRef.current.toggleFocusMode();
      setIsActive(!isActive);
    }
  };

  const startScreenShare = async () => {
    if (managerRef.current) {
      try {
        await managerRef.current.startScreenShare();
        setIsScreenSharing(true);
      } catch (error) {
        console.error('Screen sharing failed:', error);
      }
    }
  };

  return {
    isActive,
    isScreenSharing,
    toggleFocusMode,
    startScreenShare,
    manager: managerRef.current
  };
};
