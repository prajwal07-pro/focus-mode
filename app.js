const shareButton = document.getElementById('shareButton');
const videoElement = document.getElementById('videoElement');
// Get a reference to our new error message element
const errorMessage = document.getElementById('errorMessage');

shareButton.addEventListener('click', async () => {
    // Clear any previous error messages
    errorMessage.textContent = '';
    videoElement.srcObject = null;

    try {
        const stream = await navigator.mediaDevices.getDisplayMedia({
            video: true,
            audio: false
        });

        const [videoTrack] = stream.getVideoTracks();
        const settings = videoTrack.getSettings();

        if (settings.displaySurface !== 'monitor') {
            videoTrack.stop();
            // Instead of an alert, set the text content of our error element
            errorMessage.textContent = "Please share your ENTIRE SCREEN. Click the button to try again.";
            return;
        }

        videoElement.srcObject = stream;

        videoTrack.addEventListener('ended', () => {
            console.log("The user has ended the screen share.");
        });

    } catch (error) {
        // The user clicked "Cancel"
        errorMessage.textContent = "Screen share request was cancelled.";
        console.error("Error starting screen share:", error.message);
    }
});